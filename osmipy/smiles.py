import copy

import osmipy.smiles_ast
from osmipy import parser, lexer, validator, visitor
from osmipy.tokens import *


class Interpreter(visitor.NodeVisitor):
    """Visitor that gives a string representation of a smiles AST

    :param node: the node
    :type node: Chain
    """
    def __init__(self, node):
        self.node = node

    def interpret(self):
        return self.visit(self.node)

    def visit_chain(self, node):
        """

        :param node: node
        :type node: Chain
        :rtype: str
        """
        r = self.visit(node.left)
        if node.right is not None:
            if node.bond is not None:
                r += self.visit(node.bond)
            r += self.visit(node.right)

        return r

    def visit_branchedatom(self, node):
        """

        :param node: node
        :type node: BranchedAtom
        :rtype: str
        """
        r = self.visit(node.atom)
        for ringbond in node.ring_bonds:
            r += self.visit(ringbond)
        for branch in node.branches:
            r += self.visit(branch)

        return r

    def visit_atom(self, node):
        """

        :param node: node
        :type node: qcip_tools.smiles.Atom
        :rtype: str
        """
        bracketed = node.is_bracketed()
        charge = ''
        if node.charge != 0:
            charge = '+' if node.charge > 0 else '-'
            if node.charge > 1 or node.charge < -1:
                charge += str(node.charge)

        return '{}{}{}{}{}{}{}{}'.format(
            LSPAR if bracketed else '',
            node.isotope if node.isotope > 0 else '',
            node.symbol,
            node.chirality if node.chirality is not None else '',
            ('H' + (str(node.hcount) if node.hcount > 1 else '')) if node.hcount > 0 else '',
            charge,
            '{}{}'.format(COLON, node.klass) if node.klass > 0 else '',
            RSPAR if bracketed else '')

    def visit_branch(self, node):
        """

        :param node: node
        :type node: Branch
        :rtype: str
        """
        if node.bond is not None:
            r = self.visit(node.bond)
        else:
            r = ''
        r += self.visit(node.chain)

        return LPAR + r + RPAR

    def visit_ringbond(self, node):
        """

        :param node: node
        :type node: RingBond
        :rtype: str
        """
        return '{}{}{}'.format(
            '' if node.bond is None else node.bond.symbol,
            PERCENT if node.ring_id > 9 else '',
            node.ring_id)

    def visit_bond(self, node):
        """

        :param node: node
        :type node: Bond
        :rtype: str
        """
        return node.symbol if node.symbol is not None else ''


class SMILES:
    """SMILES

    :param input_: input
    :type input_: Chain|str
    """
    def __init__(self, input_='', validate=True):
        self.node = None
        if type(input_) is str:
            self.node = parser.Parser(lexer.Lexer(input_)).smiles()
        elif type(input_) is osmipy.smiles_ast.Chain:
            self.node = input_
        else:
            raise TypeError(input_)

        self.validator = validator.Validator(self.node)
        if self.node is not None and validate:
            self.validator.validate()

    def __repr__(self):
        return Interpreter(self.node).interpret()

    def add_fragment(self, other, validate=True):
        """Add another fragment at the end using a DOT bond

        :param other: other smile
        :type other: SMILES|Chain
        :param validate: (re)validate the new fragment (otherwise, ``atom_ids`` is not updated!)
        :type validate: bool
        :rtype: SMILES
        """

        ns = copy.deepcopy(self)
        c = ns.node
        while c.right is not None:
            c = c.right

        if type(other) is SMILES:
            node = copy.deepcopy(other.node)
        else:
            node = copy.deepcopy(other)

        c.right = node
        node.parent = c
        c.bond = osmipy.smiles_ast.Bond('.')

        if validate:
            ns.validator.update(node)

        return ns

    def __add__(self, other):
        return self.add_fragment(other)

    def get_atom(self, atom_id):
        """Get the atom corresponding to the given id

        :param atom_id: the id
        :type atom_id: int
        :rtype: qcip_tools.smiles.Atom
        """

        return self.validator.atom_ids[atom_id]

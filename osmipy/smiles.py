import copy

import osmipy.smiles_ast
from osmipy import smiles_parser, lexer, visitor
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


class UnicityException(Exception):
    pass


class AtomIdCheckAndUpdate(visitor.ASTVisitor):
    """Visitor that

    Keeps a dictionary of the atoms id, ``atom_ids``, and check that they are all uniques.
    In other words,  it does the job of the parser for nodes that may have been constructed from scratch.

    :param node: the node
    :type node: Chain
    """
    def __init__(self, node):
        super().__init__(node)

        self.atom_ids = {}
        self.next_atom_id = 0

    def validate(self, shift_id=0):
        """Validate the AST

        :param shift_id: shift all atom id
        :type shift_id: int
        """
        if self.node is not None:
            self.atom_ids = {}
            self.next_atom_id = 0
            self._start(shift_id=shift_id)

    def visit_atom(self, node, shift_id=0, *args, **kwargs):
        """Just update the list of ids

        :param node: node
        :type node: qcip_tools.smiles.Atom
        :param shift_id: shift all atom id
        :type shift_id: int
        """

        super().visit_atom(node, *args, **kwargs)

        if node.atom_id > -1:
            if shift_id != 0:
                node.atom_id += shift_id

            if node.atom_id not in self.atom_ids:
                self.atom_ids[node.atom_id] = node
                if self.next_atom_id <= node.atom_id:
                    self.next_atom_id = node.atom_id + 1
            else:
                raise UnicityException('two atoms share the same id: {}!'.format(node.atom_id))

    def update(self, fragment):
        """Update the validator with a new fragment

        :param fragment: the new fragment
        :type fragment: Chain
        """

        validator = AtomIdCheckAndUpdate(fragment)
        validator.validate(shift_id=self.next_atom_id)
        self.atom_ids.update(validator.atom_ids)
        self.next_atom_id = validator.next_atom_id


class SMILES:
    """SMILES object

    This object is immutable.

    :param input_: input
    :type input_: Chain|str
    """
    def __init__(self, input_=''):
        self.node = None
        if type(input_) is str:
            parser_obj = smiles_parser.Parser(lexer.Lexer(input_))
            self.node = parser_obj.smiles()
            self.atom_ids = parser_obj.atom_ids
            self.next_atom_id = parser_obj.next_atom_id
        elif type(input_) is osmipy.smiles_ast.Chain:
            self.node = input_
            validator = AtomIdCheckAndUpdate(self.node)
            validator.validate()
            self.atom_ids = validator.atom_ids
            self.next_atom_id = validator.next_atom_id
        else:
            raise TypeError(input_)

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
            updater = AtomIdCheckAndUpdate(node)
            updater.validate(shift_id=self.next_atom_id)
            ns.atom_ids.update(updater.atom_ids)

        return ns

    def __add__(self, other):
        return self.add_fragment(other)

    def get_atom(self, atom_id):
        """Get the atom corresponding to the given id

        :param atom_id: the id
        :type atom_id: int
        :rtype: qcip_tools.smiles.Atom
        """

        return self.atom_ids[atom_id]

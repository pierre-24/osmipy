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

    def __call__(self):
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
        bracketed = node.bracketed
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


class DuplicateAtomIdException(Exception):
    pass


class AtomIdCheckAndUpdate(visitor.ASTVisitor):
    """Visitor that keeps a dictionary of the atoms id, ``atom_ids``, and check that they are all uniques.
    In other words, it does the job of the parser for nodes that may have been constructed from scratch.

    If an atom id is negative, sets it to a positive and unique value.

    :param node: the node
    :type node: Chain
    """
    def __init__(self, node):
        super().__init__(node)

        self.atom_ids = {}
        self.atoms_no_id = []
        self.next_atom_id = 0

    def __call__(self, shift_id=0):
        """Does the job

        :param shift_id: shift all atom id
        :type shift_id: int
        """

        if self.node is not None:
            self._start(shift_id=shift_id)

            for i in self.atoms_no_id:
                i.atom_id = self.next_atom_id
                self.next_atom_id += 1

    def visit_atom(self, node, shift_id=0, *args, **kwargs):
        """Just update the list of ids, or mark the atom in order for it to get an id.

        :param node: node
        :type node: osmipy.smiles_ast.Atom
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
                raise DuplicateAtomIdException('two atoms share the same id: {}!'.format(node.atom_id))
        else:
            self.atoms_no_id.append(node)


class SMILES:
    """SMILES object

    This object is immutable.

    :param input_: input
    :type input_: Chain|str
    """

    #: Interpret the AST and gives a string
    interpreter = Interpreter
    #: Check (for duplicates) and update (shift or set) the atom id
    id_checker = AtomIdCheckAndUpdate

    def __init__(self, input_=''):
        self.node = None
        if type(input_) is str:
            parser_obj = smiles_parser.Parser(lexer.Lexer(input_))
            self.node = parser_obj.smiles()
            self.atom_ids = parser_obj.atom_ids
            self.next_atom_id = parser_obj.next_atom_id
        elif type(input_) is osmipy.smiles_ast.Chain:
            self.node = input_
            validator = self.id_checker(self.node)
            validator()
            self.atom_ids = validator.atom_ids
            self.next_atom_id = validator.next_atom_id
        else:
            raise TypeError(input_)

    def __repr__(self):
        return self.interpreter(self.node)()

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
            validator = self.id_checker(node)
            validator(shift_id=self.next_atom_id)
            ns.atom_ids.update(validator.atom_ids)

        return ns

    def __add__(self, other):
        return self.add_fragment(other)

    def get_atom(self, atom_id):
        """Get the atom corresponding to the given id

        :param atom_id: the id
        :type atom_id: int
        :rtype: osmipy.smiles_ast.Atom
        """

        return self.atom_ids[atom_id]

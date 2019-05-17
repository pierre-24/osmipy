from osmipy import visitor, smiles_ast
from osmipy.tokens import *


class SMILESASTVisitor(visitor.NodeVisitor):
    """Generic visitor for the SMILES AST

    :param node: node
    :type node: Chain
    """
    def __init__(self, node):
        self.node = node

    def __call__(self, *args, **kwargs):
        """Start the visit
        """
        self.visit(self.node, *args, **kwargs)

    def visit_smiles(self, node, *args, **kwargs):
        """

        :param node: node
        :type node: osmipy.smiles_ast.SMILES
        """
        self.visit(node.chain, *args, **kwargs)

    def visit_chain(self, node, *args, **kwargs):
        """

        :param node: node
        :type node: osmipy.smiles_ast.Chain
        """
        self.visit(node.branched_atom, *args, **kwargs)
        if node.bond is not None:
            self.visit(node.bond, *args, **kwargs)
        if node.next_chain is not None:
            self.visit(node.next_chain, *args, **kwargs)

    def visit_branchedatom(self, node, *args, **kwargs):
        """

        :param node: node
        :type node: osmipy.smiles_ast.BranchedAtom
        """

        self.visit(node.atom, *args, **kwargs)

        for ringbond in node.ring_bonds:
            self.visit(ringbond, *args, **kwargs)
        for branch in node.branches:
            self.visit(branch, *args, **kwargs)

    def visit_atom(self, node, *args, **kwargs):
        """

        :param node: node
        :type node: qcip_tools.smiles.Atom
        """
        pass

    def visit_bond(self, node, *args, **kwargs):
        """

        :param node: node
        :type node: qcip_tools.smiles.Bond
        """
        pass

    def visit_branch(self, node, *args, **kwargs):
        """

        :param node: node
        :type node: osmipy.smiles_ast.Branch
        """
        if node.bond is not None:
            self.visit(node.bond, *args, **kwargs)

        self.visit(node.chain, *args, **kwargs)

    def visit_ringbond(self, node, *args, **kwargs):
        """

        :param node: node
        :type node: osmipy.smiles_ast.RingBond
        """

        if node.bond is not None:
            self.visit(node.bond, *args, **kwargs)


class Interpreter(visitor.NodeVisitor):
    """Visitor that gives a string representation of a smiles AST

    :param node: the node
    :type node: Chain
    """
    def __init__(self, node):
        self.node = node

    def __call__(self):
        return self.visit(self.node)

    def visit_smiles(self, node):
        """

        :param node: node
        :type node: osmipy.smiles_ast.SMILES
        :rtype: str
        """

        return self.visit(node.chain)

    def visit_chain(self, node):
        """

        :param node: node
        :type node: osmipy.smiles_ast.Chain
        :rtype: str
        """
        r = self.visit(node.branched_atom)
        if node.next_chain is not None:
            if node.bond is not None:
                r += self.visit(node.bond)
            r += self.visit(node.next_chain)

        return r

    def visit_branchedatom(self, node):
        """

        :param node: node
        :type node: osmipy.smiles_ast.BranchedAtom
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
            node.configuration if node.configuration is not None else '',
            ('H' + (str(node.hcount) if node.hcount > 1 else '')) if node.hcount > 0 else '',
            charge,
            '{}{}'.format(COLON, node.klass) if node.klass > 0 else '',
            RSPAR if bracketed else '')

    def visit_branch(self, node):
        """

        :param node: node
        :type node: osmipy.smiles_ast.Branch
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
        :type node: osmipy.smiles_ast.RingBond
        :rtype: str
        """
        return '{}{}{}'.format(
            '' if node.bond is None else node.bond.symbol,
            PERCENT if node.ring_id > 9 else '',
            node.ring_id)

    def visit_bond(self, node):
        """

        :param node: node
        :type node: osmipy.smiles_ast.Bond
        :rtype: str
        """
        return node.symbol if node.symbol is not None else ''


class DuplicateAtomIdException(Exception):
    pass


class AtomIdCheckAndUpdate(SMILESASTVisitor):
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

    def __call__(self, shift_id=0, remove_duplicate=False):
        """Does the job

        :param shift_id: shift all atom id
        :type shift_id: int
        """

        if self.node is not None:
            super().__call__(shift_id=shift_id, remove_duplicate=remove_duplicate)

            if self.next_atom_id == 0:
                self.next_atom_id = shift_id

            for i in self.atoms_no_id:
                i.atom_id = self.next_atom_id
                self.atom_ids[self.next_atom_id] = i
                self.next_atom_id += 1

    def visit_atom(self, node, shift_id=0, remove_duplicate=False, *args, **kwargs):
        """Just update the list of ids, or mark the atom in order for it to get an id.

        :param node: node
        :type node: osmipy.smiles_ast.Atom
        :param shift_id: shift all atom id
        :type shift_id: int
        :param remove_duplicate: give a new id to duplicate instead of raising an error
        :type remove_duplicate: bool
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
                if not remove_duplicate:
                    raise DuplicateAtomIdException('two atoms share the same id: {}!'.format(node.atom_id))
                else:
                    self.atoms_no_id.append(node)
        else:
            self.atoms_no_id.append(node)


class TreeView(SMILESASTVisitor):
    """Generic visitor for the SMILES AST

    :param node: node
    :type node: Chain
    """

    NODE_CONT = '├──'
    NODE_END  = '└──'  # noqa
    NODE_TR   = '|  '  # noqa
    NODE_TRN  = '   '  # noqa

    def __init__(self, node):
        super().__init__(node)

    def __call__(self, *args, **kwargs):
        """Start the visit
        """
        return self.visit(self.node, *args, **kwargs)

    @staticmethod
    def get_levels(levels, current_is_last=False, bond=None):
        r = ''
        for i in levels:
            if i:
                r += TreeView.NODE_TR
            else:
                r += TreeView.NODE_TRN
        if current_is_last:
            r += TreeView.NODE_END
        else:
            r += TreeView.NODE_CONT

        return r

    def visit_smiles(self, node, *args, **kwargs):
        """

        :param node: node
        :type node: osmipy.smiles_ast.SMILES
        """

        r = 'SMILES'
        if node.chain:
            r += '\n' + self.visit(node.chain, levels=[])

        return r

    def visit_chain(self, node, *args, **kwargs):
        """

        :param node: node
        :type node: osmipy.smiles_ast.Chain
        """

        levels = kwargs.get('levels', [])
        r = TreeView.get_levels(levels, current_is_last=node.next_chain is None)

        if node.parent is not None:
            if type(node.parent) in [smiles_ast.Chain, smiles_ast.Branch]:
                if node.parent.bond is not None:
                    r += '[{}]'.format(str(node.parent.bond))

        levels.append(node.next_chain is not None)
        r += self.visit(node.branched_atom, levels=levels)
        levels.pop()

        if node.next_chain is not None:
            r += '\n' + self.visit(node.next_chain, levels=levels)

        return r

    def visit_branchedatom(self, node, *args, **kwargs):
        levels = kwargs.get('levels', [])
        r = str(node.atom)

        for i, n in enumerate(node.branches):
            levels.append(not(i == (len(node.branches) - 1) and len(node.ring_bonds) == 0))
            r += '\n' + self.visit(n, levels=levels)
            levels.pop()

        for i, n in enumerate(node.ring_bonds):
            levels.append(i == (len(node.ring_bonds) - 1))
            r += '\n' + self.visit(n, levels=levels)
            levels.pop()

        return r

    def visit_branch(self, node, *args, **kwargs):
        levels = kwargs.get('levels', [])

        r = TreeView.get_levels(levels[:-1], current_is_last=not levels[-1])
        r += 'BRANCH'

        if node.chain:
            r += '\n' + self.visit(node.chain, levels=levels)

        return r

    def visit_ringbond(self, node, *args, **kwargs):
        levels = kwargs.get('levels', [])
        r = TreeView.get_levels(levels[:-1], current_is_last=levels[-1]) + '>'
        r += '{}'.format(node.ring_id)
        if node.bond is not None:
            r += '[{}]'.format(node.bond)
        return r

import warnings

from osmipy.tokens import *


class AST:
    """AST element
    """
    def __init__(self):
        self.parent = None

    @property
    def parents(self):
        """Iterate over all parents
        """

        p = self.parent
        while p is not None:
            yield p
            p = p.parent

    @property
    def children(self):
        """Iterate over children
        """
        return
        yield  # yield nothing

    @property
    def descendants(self):
        """Iterate over all tag children, recursively
        """

        for c in self.children:
            yield c
            if hasattr(c, 'descendants'):
                yield from c.descendants


class Chain(AST):
    """AST element: Chain (``chain``)

    :param left: left value
    :type left: BranchedAtom
    :param right: right value
    :type right: Chain
    :param bond: bond
    :type bond: Bond
    """
    def __init__(self, left, right=None, bond=None):
        super().__init__()
        self.left = left
        self.right = right
        self.bond = bond

        self.left.parent = self

        if self.right is not None:
            self.right.parent = self
        if self.bond is not None:
            self.bond.parent = self

    @property
    def children(self):
        yield self.left
        if self.bond is not None:
            yield self.bond
        if self.right is not None:
            yield self.right


class BranchedAtom(AST):
    """AST element: BranchedAtom (``branched_atom``)

    :param atom: atom
    :type atom: qcip_tools.smiles.Atom
    :param ring_bonds: ring_bonds
    :type ring_bonds: list of RingBond
    :param branches: branches
    :type branches: list of Branch
    """

    def __init__(self, atom, ring_bonds=None, branches=None):
        super().__init__()
        self.atom = atom
        self.ring_bonds = [] if ring_bonds is None else ring_bonds
        self.branches = [] if branches is None else branches

        self.atom.parent = self

        for r in self.ring_bonds:
            r.parent = self

        for b in self.branches:
            b.parent = self

    @property
    def children(self):
        yield self.atom
        for b in self.branches:
            yield b
        for rb in self.ring_bonds:
            yield rb

    def is_organic(self):
        """

        :rtype: bool
        """
        return self.atom.is_organic()

    def is_aromatic(self):
        """

        :rtype: bool
        """
        return self.atom.is_aromatic()

    def implicit_hcount(self):
        """Determine the implicit hydrogen count by making the difference between the sum of all bond orders and
        the "normal valence" of the atom.

        **Only for the atom in the "organic" subset, not bracketed and validated!**

        :rtype: int
        """
        def bond_order(other, is_aromatic, look_left=False):
            is_other_aromatic = False
            bond = other.bond
            if type(other) is RingBond:
                if bond is None and other.target is not None:
                    bond = other.target.ring_bonds[
                        next(i for i, a in enumerate(other.target.ring_bonds) if other.parent == a.target)].bond
            if bond is None:
                if is_aromatic:
                    if type(other) is RingBond:
                        is_other_aromatic = other.target.is_aromatic()
                    elif type(other) is Branch:
                        is_other_aromatic = other.chain.left.is_aromatic()
                    elif type(other) is Chain:
                        if not look_left:
                            is_other_aromatic = other.right.left.is_aromatic()
                        else:
                            is_other_aromatic = other.left.is_aromatic()
                bo = 1.5 if (is_aromatic and is_other_aromatic) else 1
            else:
                bo = bond.bond_order()

            return bo

        if not self.is_organic() or self.atom.is_bracketed():
            return -1
        else:
            aromatic = self.is_aromatic()
            total_bonds = 0

            if self.parent.right is not None:
                total_bonds = bond_order(self.parent, aromatic)

            if self.parent.parent is not None:
                total_bonds += bond_order(self.parent.parent, aromatic, look_left=True)

            for i in self.ring_bonds:
                total_bonds += bond_order(i, aromatic)
            for i in self.branches:
                total_bonds += bond_order(i, aromatic)

            if total_bonds != int(total_bonds):
                warnings.warn('fractional bond_order of {:.1f}'.format(total_bonds), category=RuntimeWarning)

            total_bonds = int(total_bonds)

            normal_valences = NORMAL_VALENCES[self.atom.atom_symbol()]
            for v in normal_valences:
                if v < total_bonds:
                    continue
                else:
                    return v - total_bonds

            return 0

    def number_of_neighbours(self):
        """Return the number of neighbours that the atom have

        :rtype: int
        """
        n = 0

        if self.parent.right is not None:
            n += 1
        if self.parent.parent is not None:
            n += 1

        n += len(self.branches) + len(self.ring_bonds)

        if self.atom.is_bracketed():
            n += self.atom.hcount
        else:
            n += self.implicit_hcount()

        return n


class Branch(AST):
    """AST element: Branch (``branch``)

    :param chain: chain
    :type chain: Chain
    :param bond: bond
    :type bond: Bond
    """
    def __init__(self, chain, bond=None):
        super().__init__()
        self.chain = chain
        self.bond = bond

        self.chain.parent = self

        if self.bond is not None:
            self.bond.parent = self

    @property
    def children(self):
        if self.bond is not None:
            yield self.bond

        yield self.chain


class RingBond(AST):
    """AST element: RingBond (``(bond | DOT) ? ring_id``)

    :param ring_id: atom identifier (<100)
    :type ring_id: int
    :param bond: bond
    :type bond: Bond
    :param target: the other atom
    :type target: BranchedAtom
    """
    def __init__(self, ring_id, bond=None, target=None):
        super().__init__()
        self.ring_id = ring_id
        self.bond = bond
        self.target = target

        if self.bond is not None:
            self.bond.parent = self


class Atom(AST):
    """AST element: Atom (``atom``)

    :param symbol: atomic symbol
    :type symbol: str
    :param isotope: isotope
    :type isotope: int
    :param chirality: chirality symbol
    :type chirality: str
    :param hcount: hydrogen count
    :type hcount: int
    :param charge: charge
    :type charge: int
    :param klass: class
    :type klass: int
    """
    def __init__(self, symbol, isotope=0, chirality=None, hcount=0, charge=0, klass=0, atom_id=-1):
        super().__init__()
        self.symbol = symbol
        self.isotope = isotope
        self.chirality = chirality
        self.hcount = hcount
        self.charge = charge
        self.klass = klass
        self.atom_id = atom_id

    def is_bracketed(self):
        """Should this atom be bracketed?

        :rtype: bool
        """
        return self.isotope > 0 or \
            self.chirality is not None or \
            self.hcount != 0 or \
            self.charge != 0 or \
            self.klass > 0 or \
            self.symbol not in ORGANIC_SUBSET + [WILDCARD]

    def is_organic(self):
        """

        :rtype: bool
        """
        return self.symbol in ORGANIC_SUBSET

    def is_aromatic(self):
        """

        :rtype: bool
        """
        return self.symbol in AROMATIC_SYMBOLS

    def atom_symbol(self):
        if self.is_aromatic():
            return self.symbol.title()
        else:
            return self.symbol


class Bond(AST):
    """AST element: Bond (``bond | DOT``)

    :param symbol: bond symbol
    :type symbol: str
    """
    def __init__(self, symbol=None):
        super().__init__()
        self.symbol = symbol

    def bond_order(self):
        if self.symbol is None:
            return 1
        try:
            return BOND_ORDER[self.symbol]
        except KeyError:
            return 0

    def same_category(self, other):
        """Determine if two bonds are of the same category, used for ring bonds.

        The categories are:

        + If one of the bond is implicit: simple or aromatic for the other (``-``,``/``,``\\``, ``:``) ;
        + If the other is explicit:

          + Around double bond (``/``,``\\``) ;
          + Same symbol.

        :param other: the other bond
        :type other: str|Bond
        :rtype: bool
        """

        if other is None:
            symbol = None
        elif type(other) is Bond:
            symbol = other.symbol
        elif type(other) is str:
            symbol = other
        else:
            raise TypeError(other)

        if symbol is None:
            return self.symbol in ['-', '/', '\\', ':']
        if self.symbol in ['/', '\\'] and symbol in ['/', '\\']:
            return True
        else:
            return self.symbol == symbol

    @staticmethod
    def invsign(s):
        return '\\' if s == '/' else '/'

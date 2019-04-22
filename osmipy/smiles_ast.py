import warnings

from osmipy.tokens import *


class NotAChildError(ValueError):
    pass


class AST:
    """AST element
    """

    def __init__(self):
        self.parent = None

    def _set_if(self, value, attr, type_):
        if type(value) is not type_:
            raise TypeError(value)

        value.parent = self
        setattr(self, attr, value)

    def _set_if_or_none(self, value, attr, type_):
        if value is None:
            setattr(self, attr, None)
        else:
            self._set_if(value, attr, type_)

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

    @property
    def contents(self):
        """Return a list of all the children

        :rtype: list
        """
        return list(self.children)

    def _siblings_of(self, child, reverse=False):
        """Yield the next (or previous) children w.r.t. a given element

        :yield ValueError: if ``element`` is not a children
        :param child: a child
        :type child: AST
        :param reverse: yield previous children instead of next ones
        :type reverse: bool
        """

        contents = self.contents
        if reverse:
            contents.reverse()

        i = contents.index(child)

        if i < len(contents) - 1:
            yield from contents[i + 1:]  # not the last element
        else:
            return

    @property
    def next_siblings(self):
        """Yield next siblings on the same level of the tree"""
        if self.parent is not None:
            yield from self.parent._siblings_of(self)
        else:
            return

    @property
    def previous_siblings(self):
        """Yield previous siblings on the same level of the tree"""
        if self.parent is not None:
            yield from self.parent._siblings_of(self, reverse=True)
        else:
            return

    @property
    def previous_sibling(self):
        """Get next siblings on the same level of the tree
        """

        try:
            return next(self.previous_siblings)
        except StopIteration:
            return None

    @property
    def next_sibling(self):
        """Get previous siblings on the same level of the tree
        """

        try:
            return next(self.next_siblings)
        except StopIteration:
            return None

    def signal_remove(self):
        """Signal that the node will be removed, so that ring bonds get cleaned up

        The signal is broadcasted to the children (and so all)
        """

        for c in self.children:
            c.signal_remove()

    def _remove_child(self, child):
        raise NotImplementedError('_remove_child()')

    def remove(self, signal_remove=True):
        """Remove itself from the parent

        :param signal_remove: use ``signal_remove()`` before
        :type signal_remove: bool
        """

        if self.parent is None:
            raise AttributeError('`self` has no parent')

        if signal_remove:
            self.signal_remove()

        return self.parent._remove_child(self)

    def _replace_child_with(self, child, node):
        """Replace a children by another (at the same position)
        """

        raise NotImplementedError('_replace_child()')

    def replace_with(self, node, signal_remove=True):
        """Replace itself by another node

        :param node: the node with which ``self`` will be replaced
        :type node: AST
        :param signal_remove: use ``signal_remove()`` before
        :type signal_remove: bool
        """

        if self.parent is None:
            raise AttributeError('`self` has no parent')

        if signal_remove:
            self.signal_remove()

        self.parent._replace_child_with(self, node)


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

        self._left = None
        self._right = None
        self._bond = None

        # set
        self.left = left
        self.right = right
        self.bond = bond

    @property
    def left(self):
        return self._left

    @left.setter
    def left(self, value):
        self._set_if(value, '_left', BranchedAtom)

    @property
    def right(self):
        return self._right

    @right.setter
    def right(self, value):
        self._set_if_or_none(value, '_right', Chain)

    @property
    def bond(self):
        return self._bond

    @bond.setter
    def bond(self, value):
        self._set_if_or_none(value, '_bond', Bond)

    @property
    def children(self):
        yield self._left
        if self._bond is not None:
            yield self._bond
        if self._right is not None:
            yield self._right

    def _remove_child(self, child):
        if type(child) is BranchedAtom:
            raise ValueError('cannot remove a `BranchedAtom`')
        elif child == self._bond:
            self._bond = None
        elif child == self._right:
            self._right = None
        else:
            raise NotAChildError(child)

    def _replace_child_with(self, child, node):
        if child == self._left:
            self.left = node
        elif child == self._bond:
            self.bond = node
        elif child == self._right:
            self.right = node
        else:
            raise NotAChildError(child)


class NotOrganicException(Exception):
    pass


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

        self._atom = None

        # set
        self.atom = atom
        self.branches = []
        self.ring_bonds = []

        if ring_bonds is not None:
            for i in ring_bonds:
                self.append_ring_bond(i)

        if branches is not None:
            for i in branches:
                self.append_branch(i)

    @property
    def atom(self):
        return self._atom

    @atom.setter
    def atom(self, value):
        self._set_if(value, '_atom', Atom)

    def append_branch(self, branch):
        """Add a branch

        :type branch: Branch
        """

        if type(branch) is not Branch:
            raise TypeError(branch)

        branch.parent = self
        self.branches.append(branch)

    def append_ring_bond(self, ring_bond):
        """Add a ring bond

        :type ring_bond: RingBond
        """

        if type(ring_bond) is not RingBond:
            raise TypeError(ring_bond)

        ring_bond.parent = self
        self.ring_bonds.append(ring_bond)

    @property
    def children(self):
        yield self._atom
        for rb in self.ring_bonds:
            yield rb
        for b in self.branches:
            yield b

    def _remove_child(self, child):
        if type(child) is Atom:
            raise ValueError('cannot remove an `Atom`')
        elif type(child) is Branch:
            try:
                i = self.branches.index(child)
                del self.branches[i]
            except ValueError:
                raise NotAChildError(child)
        elif type(child) is RingBond:
            try:
                i = self.ring_bonds.index(child)
                del self.ring_bonds[i]  # both end will be deleted
            except ValueError:
                raise NotAChildError(child)
        else:
            raise NotAChildError(child)

    def _replace_child_with(self, child, node):
        if child == self._atom:
            self.atom = node
        elif type(child) is Branch:
            if type(node) is not Branch:
                raise TypeError(node)
            try:
                i = self.branches.index(child)
                self.branches[i] = node
                node.parent = self
            except ValueError:
                raise NotAChildError(child)
        elif type(child) is RingBond:
            if type(node) is not RingBond:
                raise TypeError(node)
            try:
                i = self.ring_bonds.index(child)
                self.ring_bonds[i] = node
                node.parent = self
            except ValueError:
                raise NotAChildError(child)
        else:
            raise NotAChildError(child)

    def _insert_child(self, child, node, after=False):
        if type(child) is Atom:
            raise ValueError('cannot insert before/after an `Atom`')
        elif type(child) is Branch:
            if type(node) is not Branch:
                raise TypeError(node)
            try:
                i = self.branches.index(child)
                if not after:
                    self.branches.insert(i, node)
                else:
                    self.branches.insert(i + 1, node)
                node.parent = self
            except ValueError:
                raise NotAChildError(child)
        elif type(child) is RingBond:
            if type(node) is not RingBond:
                raise TypeError(node)
            try:
                i = self.ring_bonds.index(child)
                if not after:
                    self.ring_bonds.insert(i, node)
                else:
                    self.ring_bonds.insert(i + 1, node)
                node.parent = self
            except ValueError:
                raise NotAChildError(child)
        else:
            raise NotAChildError(child)

    def _bond_order_with(self, neighbour, look_left=False):
        """Determine the bond order with a given neighbour

        :param neighbour: the neighbour to ``self``
        :type neighbour: Chain|RingBond|Branch
        :param look_left: in a chain, look in the ``left`` branch
        :type look_left: bool
        :rtype: int
        """

        is_other_aromatic = False
        bond = neighbour.bond

        if type(neighbour) is RingBond:
            if bond is None and neighbour.target is not None:
                bond = neighbour.target.ring_bonds[
                    next(i for i, a in enumerate(neighbour.target.ring_bonds) if neighbour.parent == a.target)].bond

        if bond is None:
            if self._atom.aromatic:
                if type(neighbour) is RingBond:
                    is_other_aromatic = neighbour.target.atom.aromatic
                elif type(neighbour) is Branch:
                    is_other_aromatic = neighbour.chain.left.atom.aromatic
                elif type(neighbour) is Chain:
                    if not look_left:
                        is_other_aromatic = neighbour.right.left.atom.aromatic
                    else:
                        is_other_aromatic = neighbour.left.aromatic
            bo = 1.5 if (self._atom.aromatic and is_other_aromatic) else 1
        else:
            bo = bond.bond_order

        return bo

    @property
    def implicit_hcount(self):
        """
        Determine the implicit hydrogen count by making the difference between the sum of all bond orders and
        the "normal valence" of the atom.

        **Only works for the atom in the "organic" subset**

        :rtype: int
        """

        if not self._atom.organic or self._atom.bracketed:
            raise NotOrganicException(self._atom)

        total_bonds = 0

        if self.parent.parent is not None:
            total_bonds = self._bond_order_with(self.parent.parent, look_left=True)

        if self.parent.right is not None:
            total_bonds += self._bond_order_with(self.parent)

        for i in self.ring_bonds:
            total_bonds += self._bond_order_with(i)
        for i in self.branches:
            total_bonds += self._bond_order_with(i)

        if total_bonds != int(total_bonds):
            warnings.warn('fractional bond_order of {:.1f}'.format(total_bonds), category=RuntimeWarning)

        total_bonds = int(total_bonds)

        normal_valences = NORMAL_VALENCES[self._atom.atom_symbol()]

        for v in normal_valences:
            if v < total_bonds:
                continue
            else:
                return v - total_bonds

        return 0  # hypervalent atom

    @property
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

        if self._atom.bracketed():
            n += self._atom.hcount
        else:
            n += self.implicit_hcount

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

        self._chain = None
        self._bond = None

        # set
        self.chain = chain
        self.bond = bond

    @property
    def bond(self):
        return self._bond

    @bond.setter
    def bond(self, value):
        self._set_if_or_none(value, '_bond', Bond)

    @property
    def chain(self):
        return self._chain

    @chain.setter
    def chain(self, value):
        self._set_if(value, '_chain', Chain)

    @property
    def children(self):
        if self._bond is not None:
            yield self._bond

        yield self._chain

    def _remove_child(self, child):
        if child == self._chain:
            self._chain = None
        elif child == self._bond:
            self._bond = None
        else:
            raise NotAChildError(child)

    def _replace_child_with(self, child, node):
        if child == self._chain:
            self.chain = node
        elif child == self._bond:
            self.bond = node
        else:
            raise NotAChildError(child)

    def _insert_child(self, child, node, after=False):
        raise NotImplementedError('_insert_child()')

    def insert_before(self, node):
        """Insert a chain before the current one

        :param node: the chain which will be inserted
        :type node: Chain
        """

        if self.parent is None:
            raise AttributeError('`self` has no parent')

        self.parent._insert_child(self, node)

    def insert_after(self, node):
        """Insert a chain after the current one

        :param node: the chain which will be inserted
        :type node: Chain
        """

        if self.parent is None:
            raise AttributeError('`self` has no parent')

        self.parent._insert_child(self, node, after=True)


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

        self._bond = None
        self._target = None

        # set
        self.ring_id = ring_id
        self.bond = bond
        self.target = target

    @property
    def bond(self):
        return self._bond

    @bond.setter
    def bond(self, value):
        self._set_if_or_none(value, '_bond', Bond)

    @property
    def target(self):
        return self._target

    @target.setter
    def target(self, value):
        if value is None:
            self._target = None
        else:
            if type(value) is not BranchedAtom:
                raise TypeError(value)

            self._target = value

    @property
    def children(self):
        if self._bond is not None:
            yield self._bond

        return

    def signal_remove(self):
        """Just remove the other ring bond in case of deletion
        """

        try:
            i = next(i for i, a in enumerate(self._target.ring_bonds) if a._target == self.parent)
            del self._target.ring_bonds[i]
        except:  # something was wrong (parent is already gone, or ring bond was already deleted) ...
            pass  # ... Never mind.

        super().signal_remove()

    def _remove_child(self, child):
        if child == self._bond:
            self._bond = None
        else:
            raise NotAChildError(child)

    def _replace_child_with(self, child, node):
        if child == self._bond:
            self.bond = node
        else:
            raise NotAChildError(child)

    def _insert_child(self, child, node, after=False):
        raise NotImplementedError('_insert_child()')

    def insert_before(self, node):
        """Insert a ring bond before the current one

        :param node: the ring bond which will be inserted
        :type node: RingBond
        """

        if self.parent is None:
            raise AttributeError('`self` has no parent')

        self.parent._insert_child(self, node)

    def insert_after(self, node):
        """Insert a ring bond after the current one

        :param node: the ring bond which will be inserted
        :type node: RingBond
        """

        if self.parent is None:
            raise AttributeError('`self` has no parent')

        self.parent._insert_child(self, node, after=True)


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

        # set
        self.symbol = symbol
        self.isotope = isotope
        self.chirality = chirality
        self.hcount = hcount
        self.charge = charge
        self.klass = klass
        self.atom_id = atom_id

    @property
    def contents(self):
        raise AttributeError('`{}` object has no `contents`'.format(type(self).__name__))

    @property
    def bracketed(self):
        """Should this atom be bracketed?

        :rtype: bool
        """
        return self.isotope > 0 or \
            self.chirality is not None or \
            self.hcount != 0 or \
            self.charge != 0 or \
            self.klass > 0 or \
            self.symbol not in ORGANIC_SUBSET + [WILDCARD]

    @property
    def organic(self):
        """

        :rtype: bool
        """
        return self.symbol in ORGANIC_SUBSET

    @property
    def aromatic(self):
        """

        :rtype: bool
        """
        return self.symbol in AROMATIC_SYMBOLS

    def atom_symbol(self):
        if self.aromatic:
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

        # set
        self.symbol = symbol

    @property
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

    @property
    def contents(self):
        raise AttributeError('`{}` object has no `contents`'.format(type(self).__name__))

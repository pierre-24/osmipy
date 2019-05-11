from osmipy import lexer
from osmipy.smiles_ast import Chain, BranchedAtom, Branch, RingBond, Atom, Bond, SMILES
from osmipy.tokens import *


class ParserException(Exception):
    def __init__(self, token, msg):
        super().__init__('parser error at position {} [{}]: {}'.format(token.position, repr(token), msg))
        self.token = token
        self.message = msg


class Parser:
    """Parser (generate and AST from the tokens).

    :type lexer: Lexer
    :param lexer: The lexer
    """

    def __init__(self, lexer):
        self.lexer = lexer
        self.tokenizer = lexer.tokenize()
        self.current_token = None
        self.previous_tokens = []
        self.use_previous = 0

        self.next_atom_id = 0
        self.atom_ids = {}

        self._ring_ids = {}
        self._ring_pairs_pid = []
        self.ring_bond_pairs = []

        self.next()

    def eat(self, token_type):
        """Consume the token if of the right type

        :param token_type: the token type
        :type token_type: str
        :raise ParserException: if not of the correct type
        """
        if self.current_token.type == token_type:
            self.next()
        else:
            raise ParserException(self.current_token, 'token must be {}'.format(token_type))

    def next(self):
        """Get the next token
        """

        try:
            if self.current_token is not None:
                self.previous_tokens.append(self.current_token)
            self.current_token = next(self.tokenizer)
        except StopIteration:
            self.current_token = Token(EOF, None)

    def bracket_atom(self):
        """

        :rtype: qcip_tools.smiles.Atom
        """

        self.eat(LSPAR)

        isotope = 0
        chirality = None
        hcount = 0
        charge = 0
        klass = -1

        # isotope
        while self.current_token.type == DIGIT:
            isotope = isotope * 10 + self.current_token.value
            self.next()

        # symbol
        if self.current_token.type != LETTER:
            raise ParserException(self.current_token, 'expected ATOM in bracket_atom')

        symbol = self.atomic_symbol()

        # chirality
        if self.current_token.type == AT:
            chirality = '@'
            self.next()
            # TODO: extend chirality at that point!
            if self.current_token.type == AT:
                chirality += '@'
                self.next()

        # hcount
        if self.current_token.type == LETTER:
            if self.current_token.value != 'H':
                raise ParserException(self.current_token, 'expected hydrogen in hcount')
            self.next()
            hcount = 1

            if self.current_token.type == DIGIT:
                hcount = self.current_token.value
                self.next()

        # charge
        if self.current_token.type in [PLUS, MINUS]:
            sign = -1 if self.current_token.type == MINUS else 1
            self.next()
            charge = 0
            i = 0
            while self.current_token.type == DIGIT and i < 2:
                charge = charge * 10 + self.current_token.value
                i += 1
                self.next()

            if charge == 0:
                charge = 1

            charge *= sign

        # class
        if self.current_token.type == COLON:
            self.next()
            if self.current_token.type != DIGIT:
                raise ParserException(self.current_token, 'expected digit for class')

            klass = 0

            while self.current_token.type == DIGIT:
                klass = klass * 10 + self.current_token.value
                self.next()

        self.eat(RSPAR)

        return Atom(symbol=symbol, isotope=isotope, chirality=chirality, hcount=hcount, charge=charge, klass=klass)

    def atomic_symbol(self):
        symbol = self.current_token.value
        self.next()

        if self.current_token.type == LETTER:
            n_symbol = symbol + self.current_token.value
            if n_symbol in TOT_SYMBOLS:
                self.next()
                return n_symbol

        if symbol not in TOT_SYMBOLS:
            raise ParserException(self.current_token, '{} is not a valid symbol'.format(symbol))

        return symbol

    def atom(self):
        """

        :rtype: qcip_tools.smiles.Atom
        """

        if self.current_token.type == LSPAR:
            atom = self.bracket_atom()
        elif self.current_token.type == LETTER:
            if self.current_token.value not in ORGANIC_SUBSET:
                raise ParserException(self.current_token.position, '{} should be bracketed!')
            atom = Atom(symbol=self.atomic_symbol())
        elif self.current_token.type == WILDCARD:
            atom = Atom(symbol=self.current_token.value)
            self.next()
        else:
            raise ParserException(self.current_token, 'unexpected token in atom')

        atom.atom_id = self.next_atom_id
        self.atom_ids[self.next_atom_id] = atom
        self.next_atom_id += 1

        return atom

    def branch(self):
        """

        :rtype: osmipy.smiles_ast.Branch
        """

        self.eat(LPAR)

        bond = None

        if self.current_token.type in BONDS_TYPE + [DOT]:
            bond = Bond(self.current_token.value)
            self.next()

        chain = self.chain()

        self.eat(RPAR)

        return Branch(chain=chain, bond=bond)

    def ring_bond(self):
        """

        :rtype: osmipy.smiles_ast.RingBond
        """

        ring_id = 0

        if self.current_token.type == PERCENT:  # > 10
            token_percent = self.current_token
            self.next()

            if self.current_token.type != DIGIT:
                raise ParserException(token_percent, 'expected DIGIT in ringbond')

            i = 0
            while self.current_token.type == DIGIT and i < 2:
                ring_id = ring_id * 10 + self.current_token.value
                i += 1
                self.next()
        elif self.current_token.type == DIGIT:  # < 10
            ring_id = self.current_token.value
            self.next()
        else:
            raise ParserException(self.current_token, 'expected PERCENT or DIGIT in ring_id')

        node = RingBond(ring_id=ring_id)
        return node

    def _consolidate_ring_bonds(self, ring_bonds):
        """Connect the two ends of a ring bond.

        The first time the ring id is met, it stores the ring bond into ``self.ring_ids``.
        The second times, it connects the two ends, but check first:

        + That the bonds are of the same type (ex: ``C(-1)CC=1`` is not allowed, but note that ``C1C=1`` is) ;
        + That it does not bound to the same atom (ex: ``C11`` is not allowed) ;
        + That the same bond is not defined twice (ex: ``C12CC12`` is not allowed) ;
        + That it does not bond to a direct pair (ex: ``C1C1`` is not allowed).

        :param ring_bonds: list of ring bonds to consolidate
        :type ring_bonds: list[osmipy.smiles_ast.RingBond]
        """

        for rb in ring_bonds:
            i = rb.ring_id
            if i not in self._ring_ids:
                self._ring_ids[rb.ring_id] = rb  # store 'til next time
            else:  # connect
                other_rb = self._ring_ids[i]
                if rb.bond is not None and other_rb.bond is not None:
                    if not rb.bond.same_category(other_rb.bond):
                        raise ParserException(
                            self.current_token, 'ring id {}: not the same type of bond at the two ends'.format(i))

                if id(rb.parent) == id(other_rb.parent):
                    raise ParserException(self.current_token, 'ring id {}: bond to same atom'.format(i))

                pair = (id(rb.parent), id(other_rb.parent))
                if pair[0] > pair[1]:
                    pair = tuple(reversed(pair))

                if pair in self._ring_pairs_pid:
                    raise ParserException(self.current_token, 'ring id {}: this bond is already defined'.format(i))
                else:
                    self._ring_pairs_pid.append(pair)

                # if everything is ok, do the connection
                rb.target, other_rb.target = other_rb.parent, rb.parent

                del self._ring_ids[i]
                self.ring_bond_pairs.append((other_rb, rb))  # 'til the end

    def branched_atom(self):
        """
        **Only** matches

        .. code-block:: text

            atom ring_id*

        to avoid the need of backtracking.

        :rtype: osmipy.smiles_ast.BranchedAtom
        """

        atom = self.atom()
        ringbonds = []

        while self.current_token.type in [DIGIT, PERCENT]:
            ringbonds.append(self.ring_bond())

        ba = BranchedAtom(atom=atom, ring_bonds=ringbonds)

        self._consolidate_ring_bonds(ringbonds)  # ... only after setting their parents !
        return ba

    def chain(self):
        """
        Actually matches

        .. code-block:: text

            branched_atom ring_bond* branch* ((bond | DOT)? chain)?

        to avoid the need of backtracking due to bond before ``ring_id``.
        ``Branch`` and ``RingBond`` are added to ``BranchedAtom``, thus fulfilling the grammar rules.

        :rtype: osmipy.smiles_ast.Chain
        """

        left = self.branched_atom()
        bond = None

        if self.current_token.type in BONDS_TYPE + [DOT]:
            bond = Bond(self.current_token.value)
            self.next()

        right = None
        ring_bonds_to_consolidate = []

        while self.current_token.type in [PERCENT, DIGIT]:
            if bond is None and bond.symbol != DOT:
                raise ParserException(self.current_token, 'ring_id requires a bond in this position')

            ring_bond = self.ring_bond()
            ring_bond.bond = bond
            ring_bonds_to_consolidate.append(ring_bond)
            left.append_ring_bond(ring_bond)

            bond = None

            # needs to get an eventual new bond
            if self.current_token.type in BONDS_TYPE + [DOT]:
                bond = Bond(self.current_token.value)
                self.next()

        self._consolidate_ring_bonds(ring_bonds_to_consolidate)

        if bond is None:
            while self.current_token.type == LPAR:
                left.append_branch(self.branch())

            # needs to get an eventual new bond
            if self.current_token.type in BONDS_TYPE + [DOT]:
                bond = Bond(self.current_token.value)
                self.next()

        if self.current_token.type in [LETTER, LSPAR]:
            right = self.chain()
        elif bond is not None:
            raise ParserException(self.current_token, 'bond but no chain')

        return Chain(left=left, bond=bond, right=right)

    def smiles(self):
        """

        :rtype: osmipy.smiles_ast.SMILES
        """

        node = None

        if self.current_token.type != EOF:
            node = self.chain()

        self.eat(EOF)

        # check for unmatched ring bonds
        if len(self._ring_ids) != 0:
            raise ParserException(
                self.current_token,
                'unmatched ring ids left: {}'.format(','.join(str(i) for i in self._ring_ids.keys())))

        # check for direct pair
        for rb1, rb2 in self.ring_bond_pairs:
            if id(rb1.parent.parent) == id(rb2.parent.parent.parent):
                raise ParserException(self.current_token, 'ring id {}: direct pair is not allowed'.format(rb1.ring_id))

        return SMILES(node, self.atom_ids)


def parse(s):
    """Parse a string to an AST

    :param s: string to parse
    :type s: str
    :rtype: osmipy.smiles_ast.SMILES
    """

    return Parser(lexer.Lexer(s)).smiles()

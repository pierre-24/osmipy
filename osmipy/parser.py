import osmipy.tokens
from osmipy.smiles_ast import Chain, BranchedAtom, Branch, RingBond, Atom, Bond
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
        self.atom_id = 0
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
            self.current_token = osmipy.tokens.Token(EOF, None)

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
        if self.current_token.type != ATOM:
            raise ParserException(self.current_token, 'expected ATOM in bracket_atom')

        symbol = self.current_token.value
        self.next()

        # chirality
        if self.current_token.type == AT:
            chirality = '@'
            self.next()
            # TODO: extend chirality at that point!
            if self.current_token.type == AT:
                chirality += '@'
                self.next()

        # hcount
        if self.current_token.type == ATOM:
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

    def atom(self):
        """

        :rtype: qcip_tools.smiles.Atom
        """

        if self.current_token.type == LSPAR:
            atom = self.bracket_atom()
        elif self.current_token.type == ATOM:
            if self.current_token.value not in ORGANIC_SUBSET:
                raise ParserException(self.current_token.position, '{} should be bracketed!')
            atom = Atom(symbol=self.current_token.value)
            self.next()
        elif self.current_token.type == WILDCARD:
            atom = Atom(symbol=self.current_token.value)
            self.next()
        else:
            raise ParserException(self.current_token, 'unexpected token in atom')

        atom.atom_id = self.atom_id
        self.atom_id += 1

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

    def ring_id(self):
        """

        :rtype: osmipy.smiles_ast.RingBond
        """

        if self.current_token.type == PERCENT:
            token_percent = self.current_token
            self.next()

            if self.current_token.type != DIGIT:
                raise ParserException(token_percent, 'expected DIGIT in ringbond')

            number = 0
            i = 0
            while self.current_token.type == DIGIT and i < 2:
                number = number * 10 + self.current_token.value
                i += 1
                self.next()
            return RingBond(ring_id=number)
        elif self.current_token.type == DIGIT:
            node = RingBond(ring_id=self.current_token.value)
            self.next()
            return node
        else:
            raise ParserException(self.current_token, 'expected PERCENT or DIGIT in ring_id')

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
            ringbonds.append(self.ring_id())

        return BranchedAtom(atom=atom, ring_bonds=ringbonds, branches=[])

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

        while self.current_token.type in [PERCENT, DIGIT]:
            if bond is None and bond.symbol != DOT:
                raise ParserException(self.current_token, 'ring_id requires a bond in this position')

            ring_id = self.ring_id()
            ring_id.bond = bond
            ring_id.parent = left
            left.ring_bonds.append(ring_id)
            bond = None

            # needs to get an eventual new bond
            if self.current_token.type in BONDS_TYPE + [DOT]:
                bond = Bond(self.current_token.value)
                self.next()

        if bond is None:
            while self.current_token.type == LPAR:
                left.branches.append(self.branch())

            # needs to get an eventual new bond
            if self.current_token.type in BONDS_TYPE + [DOT]:
                bond = Bond(self.current_token.value)
                self.next()

        if self.current_token.type in [ATOM, LSPAR]:
            right = self.chain()
        elif bond is not None:
            raise ParserException(self.current_token, 'bond but no chain')

        return Chain(left=left, bond=bond, right=right)

    def smiles(self):
        """

        :rtype: Chain
        """

        node = None

        if self.current_token.type != EOF:
            node = self.chain()

        self.eat(EOF)

        return node

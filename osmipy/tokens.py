ATOM, BOND, DIGIT, LPAR, RPAR, LSPAR, RSPAR, PLUS, MINUS, DOT, WILDCARD, PERCENT, AT, COLON, EOF = (
    'ATOM', 'BOND', 'DIGIT', '(', ')', '[', ']', '+', '-', '.', '*', '%', '@', ':', 'EOF'
)

SYMBOLS_TR = {
    '=': BOND,
    '#': BOND,
    '$': BOND,
    '/': BOND,
    '\\': BOND,
    '(': LPAR,
    ')': RPAR,
    '[': LSPAR,
    ']': RSPAR,
    '+': PLUS,
    '-': MINUS,
    '.': DOT,
    '*': WILDCARD,
    '%': PERCENT,
    '@': AT,
    ':': COLON,
    # chain terminators:
    ' ': EOF,
    '\t': EOF,
    '\n': EOF,
    '\r': EOF,
    '\0': EOF,
}

BONDS_TYPE = [BOND, MINUS, COLON]

BOND_ORDER = {
    '.': 0,
    '-': 1,
    '=': 2,
    '#': 3,
    '$': 4,
    '/': 1,
    '\\': 1
}

ELEMENT_SYMBOLS = [
    'H', 'He',
    'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
    'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
    'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
    'Cs', 'Ba', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
    'Fr', 'Ra', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Fl', 'Lv',  # TODO: next symbols
    'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
    'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr'
]

AROMATIC_SYMBOLS = ['b', 'c', 'n', 'o', 'p', 's', 'se', 'as']

ALIPHATIC_SYMBOLS = ['B', 'C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br']

TOT_SYMBOLS = ELEMENT_SYMBOLS + AROMATIC_SYMBOLS

ORGANIC_SUBSET = AROMATIC_SYMBOLS + ALIPHATIC_SYMBOLS

NORMAL_VALENCES = {
    'B': (3,),
    'C': (4,),
    'N': (3, 5),
    'O': (2,),
    'P': (3, 5),
    'S': (2, 4, 6),
    'F': (1,),
    'Cl': (1,),
    'Br': (1,),
    'I': (1,),
    'Se': (2, 4, 6),
    'As': (3, 5)
}


class Token:
    """Token class"""
    def __init__(self, type_, value, position=-1):
        self.type = type_
        self.value = value
        self.position = position

    def __repr__(self):
        return 'Token({}, {}{})'.format(
            self.type, repr(self.value), ', {}'.format(self.position) if self.position > -1 else '')

from tests import OSmiPyTestCase

from osmipy import smiles, lexer, parser, validator, smiles_ast
from osmipy.tokens import Token as T


class SMILESTestCase(OSmiPyTestCase):

    def setUp(self):
        pass

    def test_lexer(self):
        """Test that the lexer gives the good tokens"""

        tests = [
            ('cccccc', [T(smiles.ATOM, 'c')] * 6),
            ('C1=CC=CC=C1', [
                T(smiles.ATOM, 'C'),
                T(smiles.DIGIT, 1),
                T(smiles.BOND, '='),
                T(smiles.ATOM, 'C'),
                T(smiles.ATOM, 'C'),
                T(smiles.BOND, '='),
                T(smiles.ATOM, 'C'),
                T(smiles.ATOM, 'C'),
                T(smiles.BOND, '='),
                T(smiles.ATOM, 'C'),
                T(smiles.DIGIT, 1),
            ]),
            ('N[C@](Br)(O)C', [
                T(smiles.ATOM, 'N'),
                T(smiles.LSPAR, '['),
                T(smiles.ATOM, 'C'),
                T(smiles.AT, '@'),
                T(smiles.RSPAR, ']'),
                T(smiles.LPAR, '('),
                T(smiles.ATOM, 'Br'),
                T(smiles.RPAR, ')'),
                T(smiles.LPAR, '('),
                T(smiles.ATOM, 'O'),
                T(smiles.RPAR, ')'),
                T(smiles.ATOM, 'C'),
            ])
        ]

        for s, r in tests:
            smiles_lexer = lexer.Lexer(s)
            lexed = [i for i in smiles_lexer.tokenize()]
            self.assertEqual(len(lexed) - 1, len(r))
            for i, token in enumerate(r):
                self.assertEqual(token.type, lexed[i].type, msg='{} of {}'.format(i, s))
                self.assertEqual(token.value, lexed[i].value, msg='{} of {}'.format(i, s))

            self.assertEqual(lexed[-1].type, smiles.EOF)

    def test_parser(self):
        """Test the parser (on a rather simple string)"""
        smiles_parser = parser.Parser(lexer.Lexer('C(\\Cl)=C/O'))

        # output of the parser
        r = smiles_ast.Chain(
            left=smiles_ast.BranchedAtom(
                atom=smiles_ast.Atom(symbol='C'),
                branches=[
                    smiles_ast.Branch(
                        chain=smiles_ast.Chain(left=smiles_ast.BranchedAtom(atom=smiles_ast.Atom(symbol='Cl'))),
                        bond=smiles_ast.Bond('\\'))
                ]
            ),
            right=smiles_ast.Chain(
                left=smiles_ast.BranchedAtom(atom=smiles_ast.Atom(symbol='C')),
                right=smiles_ast.Chain(left=smiles_ast.BranchedAtom(atom=smiles_ast.Atom(symbol='O'))),
                bond=smiles_ast.Bond('/')
            ),
            bond=smiles_ast.Bond('='))

        rx = smiles_parser.smiles()

        # types
        self.assertEqual(type(r), type(rx))
        self.assertEqual(type(r.left), type(rx.left))
        self.assertEqual(type(r.right), type(rx.right))
        self.assertEqual(type(r.bond), type(rx.bond))
        self.assertEqual(type(r.left.branches[0]), type(rx.left.branches[0]))
        self.assertEqual(type(r.left.branches[0].chain), type(rx.left.branches[0].chain))
        self.assertEqual(type(r.left.branches[0].bond), type(rx.left.branches[0].bond))
        self.assertEqual(type(r.left.branches[0].chain.left), type(rx.left.branches[0].chain.left))
        self.assertEqual(type(r.left.branches[0].chain.left.atom), type(rx.left.branches[0].chain.left.atom))
        self.assertEqual(type(r.left.atom), type(rx.left.atom))
        self.assertEqual(type(r.right.left), type(rx.right.left))
        self.assertEqual(type(r.right.right), type(rx.right.right))
        self.assertEqual(type(r.right.bond), type(rx.right.bond))
        self.assertEqual(type(r.right.left.atom), type(rx.right.left.atom))
        self.assertEqual(type(r.right.right.left.atom), type(rx.right.right.left.atom))

        # symbols and hcounts
        self.assertEqual(r.left.branches[0].bond.symbol, rx.left.branches[0].bond.symbol)
        self.assertEqual(r.left.branches[0].chain.left.atom.symbol, rx.left.branches[0].chain.left.atom.symbol)  # -Cl
        self.assertEqual(r.left.branches[0].chain.left.implicit_hcount(), 0)
        self.assertEqual(r.bond.symbol, rx.bond.symbol)
        self.assertEqual(r.left.atom.symbol, rx.left.atom.symbol)  # sp² C
        self.assertEqual(r.left.implicit_hcount(), 1)
        self.assertEqual(r.right.bond.symbol, rx.right.bond.symbol)
        self.assertEqual(r.right.left.atom.symbol, rx.right.left.atom.symbol)  # sp² C
        self.assertEqual(r.right.right.left.atom.symbol, rx.right.right.left.atom.symbol)  # -OH
        self.assertEqual(r.right.right.left.implicit_hcount(), 1)

        # parent-child
        self.assertEqual(r.parent, None)
        self.assertEqual(r.left.parent, r)
        self.assertEqual(r.right.parent, r)
        self.assertEqual(r.bond.parent, r)
        self.assertEqual(r.left.atom.parent, r.left)
        self.assertEqual(r.left.branches[0].chain.parent, r.left.branches[0])
        self.assertEqual(r.left.branches[0].bond.parent, r.left.branches[0])

    def test_validator(self):
        """Test the validator"""

        r = smiles_ast.Chain(
            left=smiles_ast.BranchedAtom(
                atom=smiles_ast.Atom(symbol='C'),
                ring_bonds=[smiles_ast.RingBond(ring_id=1)]),
            right=smiles_ast.Chain(
                left=smiles_ast.BranchedAtom(atom=smiles_ast.Atom(symbol='C')),
                right=smiles_ast.Chain(
                    left=smiles_ast.BranchedAtom(
                        atom=smiles_ast.Atom(symbol='C'),
                        ring_bonds=[smiles_ast.RingBond(ring_id=1)])
                )))

        rx = parser.Parser(lexer.Lexer('C1CC1')).smiles()
        validator.Validator(rx).validate()

        # types
        self.assertEqual(type(r), type(rx))
        self.assertEqual(type(r.left), type(rx.left))
        self.assertEqual(type(r.left.atom), type(rx.left.atom))
        self.assertEqual(type(r.right), type(rx.right))
        self.assertEqual(type(r.right.left), type(rx.right.left))
        self.assertEqual(type(r.right.left.atom), type(rx.right.left.atom))
        self.assertEqual(type(r.right.right.left), type(rx.right.right.left))
        self.assertEqual(type(r.right.right.left.atom), type(rx.right.right.left.atom))

        # symbols and hcount
        self.assertEqual(r.left.atom.symbol, rx.left.atom.symbol)
        self.assertEqual(rx.left.implicit_hcount(), 2)  # sp³ C
        self.assertEqual(r.right.left.atom.symbol, rx.right.left.atom.symbol)
        self.assertEqual(rx.right.left.implicit_hcount(), 2)  # sp³ C
        self.assertEqual(r.right.right.left.atom.symbol, rx.right.right.left.atom.symbol)
        self.assertEqual(rx.right.right.left.implicit_hcount(), 2)  # sp³ C

        # ring bonds
        self.assertEqual(type(r.left.ring_bonds[0]), type(rx.left.ring_bonds[0]))
        self.assertEqual(type(r.right.right.left.ring_bonds[0]), type(rx.right.right.left.ring_bonds[0]))

        self.assertEqual(rx.left.ring_bonds[0].ring_id, rx.right.right.left.ring_bonds[0].ring_id)
        self.assertEqual(rx.left.ring_bonds[0].target, rx.right.right.left)
        self.assertEqual(rx.left, rx.right.right.left.ring_bonds[0].target)

        # some wrong smiles
        wrong_smiles = [
            'C11',
            'C12CCCCC12',
            'C21CCCCC12',
            'C12C2CCC1',
            'C-1CCCCC=1'
        ]
        for s in wrong_smiles:
            node = parser.Parser(lexer.Lexer(s)).smiles()
            with self.assertRaises(validator.ValidationException, msg=s):
                validator.Validator(node).validate()

    def test_smile(self):
        """Test that everything works correctly by assuming that input and output should be the same"""
        test_smiles = [
            '[Cu+2]',
            '[2H+]',
            '[CH4:2]',
            'Oc1c(*)cccc1',
            'c1ccccc1',
            'N1CC2CCCC2CC1',
            'N[C@](Br)(O)C',
            'CC(C)C(=O)C(C)C',
            'C=0CCCCC=0',
            'C12(CCCCC1)CCCCC2',
            'Oc1ccccc1.NCCO',
            'C(/F)=C/F',
            '[NH4+].[NH4+].[O-]S(=O)(=O)[S-]',
            'c1c2c3c4cc1.Br2.Cl3.Cl4'
        ]

        for smi in test_smiles:
            smile = smiles.SMILES(smi)
            self.assertEqual(smi, repr(smile))

        # test add_fragment
        a = smiles.SMILES(test_smiles[4])
        b = smiles.SMILES(test_smiles[5])
        smile = a.add_fragment(b)
        self.assertEqual(repr(smile), repr(a) + '.' + repr(b))

        self.assertEqual(smile.get_atom(a.validator.next_id).symbol, 'N')  # ok, atoms_id updated
        self.assertEqual(len(smile.validator.ring_ids[1]), len(a.validator.ring_ids[1]) + len(b.validator.ring_ids[1]))
        self.assertEqual(len(smile.validator.ring_ids[2]), len(b.validator.ring_ids[2]))

    def test_hcount(self):
        to_test = [
            ('C', 4),
            ('c1ccccc1', 1),
            ('n1ccccc1', 0),  # but that does not works for 'n1ccccc1', which should gives 1 (see below)
            ('N1C=CC=C1', 1),  # .. In the kekule form, that works
            ('N1C=CC=CC=1', 0),
            ('N(=O)O', 0),
            ('N(=O)(=O)O', 0),
            ('S(=O)=O', 0),
            ('SCCO', 1),
            ('C(=O)O', 1)
        ]

        for smi, hcount in to_test:
            smile = smiles.SMILES(smi)
            self.assertEqual(hcount, smile.node.left.implicit_hcount(), msg=smi)

    def test_get_atom(self):
        s = smiles.SMILES('CCO')
        self.assertEqual(s.get_atom(0), s.node.left.atom)
        self.assertEqual(s.get_atom(1), s.node.right.left.atom)
        self.assertEqual(s.get_atom(2), s.node.right.right.left.atom)

from tests import OSmiPyTestCase

from osmipy import lexer, smiles_parser, smiles_ast


class ParserTestCase(OSmiPyTestCase):

    def test_parser(self):
        """Test the parser (on a rather simple string)"""

        smiles_parser_obj = smiles_parser.Parser(lexer.Lexer('C(\\Cl)=C/O'))

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

        rx = smiles_parser_obj.smiles().chain

        # types
        self.assertEqual(type(r), type(rx))
        self.assertEqual(type(r.branched_atom), type(rx.branched_atom))
        self.assertEqual(type(r.next_chain), type(rx.next_chain))
        self.assertEqual(type(r.bond), type(rx.bond))
        self.assertEqual(type(r.branched_atom.branches[0]), type(rx.branched_atom.branches[0]))
        self.assertEqual(type(r.branched_atom.branches[0].chain), type(rx.branched_atom.branches[0].chain))
        self.assertEqual(type(r.branched_atom.branches[0].bond), type(rx.branched_atom.branches[0].bond))
        self.assertEqual(type(r.branched_atom.branches[0].chain.branched_atom),
                         type(rx.branched_atom.branches[0].chain.branched_atom))
        self.assertEqual(type(r.branched_atom.branches[0].chain.branched_atom.atom),
                         type(rx.branched_atom.branches[0].chain.branched_atom.atom))
        self.assertEqual(type(r.branched_atom.atom), type(rx.branched_atom.atom))
        self.assertEqual(type(r.next_chain.branched_atom), type(rx.next_chain.branched_atom))
        self.assertEqual(type(r.next_chain.next_chain), type(rx.next_chain.next_chain))
        self.assertEqual(type(r.next_chain.bond), type(rx.next_chain.bond))
        self.assertEqual(type(r.next_chain.branched_atom.atom), type(rx.next_chain.branched_atom.atom))
        self.assertEqual(type(r.next_chain.next_chain.branched_atom.atom),
                         type(rx.next_chain.next_chain.branched_atom.atom))

        # symbols and hcounts
        self.assertEqual(r.branched_atom.branches[0].bond.symbol, rx.branched_atom.branches[0].bond.symbol)
        self.assertEqual(r.branched_atom.branches[0].chain.branched_atom.atom.symbol,
                         rx.branched_atom.branches[0].chain.branched_atom.atom.symbol)  # -Cl
        self.assertEqual(r.branched_atom.branches[0].chain.branched_atom.implicit_hcount, 0)
        self.assertEqual(r.bond.symbol, rx.bond.symbol)
        self.assertEqual(r.branched_atom.atom.symbol, rx.branched_atom.atom.symbol)  # sp² C
        self.assertEqual(r.branched_atom.implicit_hcount, 1)
        self.assertEqual(r.next_chain.bond.symbol, rx.next_chain.bond.symbol)
        self.assertEqual(r.next_chain.branched_atom.atom.symbol, rx.next_chain.branched_atom.atom.symbol)  # sp² C
        self.assertEqual(r.next_chain.next_chain.branched_atom.atom.symbol,
                         rx.next_chain.next_chain.branched_atom.atom.symbol)  # -OH
        self.assertEqual(r.next_chain.next_chain.branched_atom.implicit_hcount, 1)

        # parent-child
        self.assertEqual(r.parent, None)
        self.assertEqual(r.branched_atom.parent, r)
        self.assertEqual(r.next_chain.parent, r)
        self.assertEqual(r.bond.parent, r)
        self.assertEqual(r.branched_atom.atom.parent, r.branched_atom)
        self.assertEqual(r.branched_atom.branches[0].chain.parent, r.branched_atom.branches[0])
        self.assertEqual(r.branched_atom.branches[0].bond.parent, r.branched_atom.branches[0])

        # An example with ring bonds
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

        rx = smiles_parser.Parser(lexer.Lexer('C1CC1')).smiles().chain

        # types
        self.assertEqual(type(r), type(rx))
        self.assertEqual(type(r.branched_atom), type(rx.branched_atom))
        self.assertEqual(type(r.branched_atom.atom), type(rx.branched_atom.atom))
        self.assertEqual(type(r.next_chain), type(rx.next_chain))
        self.assertEqual(type(r.next_chain.branched_atom), type(rx.next_chain.branched_atom))
        self.assertEqual(type(r.next_chain.branched_atom.atom), type(rx.next_chain.branched_atom.atom))
        self.assertEqual(type(r.next_chain.next_chain.branched_atom), type(rx.next_chain.next_chain.branched_atom))
        self.assertEqual(type(r.next_chain.next_chain.branched_atom.atom),
                         type(rx.next_chain.next_chain.branched_atom.atom))

        # symbols and hcount
        self.assertEqual(r.branched_atom.atom.symbol, rx.branched_atom.atom.symbol)
        self.assertEqual(rx.branched_atom.implicit_hcount, 2)  # sp³ C
        self.assertEqual(r.next_chain.branched_atom.atom.symbol, rx.next_chain.branched_atom.atom.symbol)
        self.assertEqual(rx.next_chain.branched_atom.implicit_hcount, 2)  # sp³ C
        self.assertEqual(r.next_chain.next_chain.branched_atom.atom.symbol,
                         rx.next_chain.next_chain.branched_atom.atom.symbol)
        self.assertEqual(rx.next_chain.next_chain.branched_atom.implicit_hcount, 2)  # sp³ C

        # ring bonds
        self.assertEqual(type(r.branched_atom.ring_bonds[0]), type(rx.branched_atom.ring_bonds[0]))
        self.assertEqual(type(r.next_chain.next_chain.branched_atom.ring_bonds[0]),
                         type(rx.next_chain.next_chain.branched_atom.ring_bonds[0]))

        self.assertEqual(rx.branched_atom.ring_bonds[0].ring_id,
                         rx.next_chain.next_chain.branched_atom.ring_bonds[0].ring_id)
        self.assertEqual(rx.branched_atom.ring_bonds[0].target, rx.next_chain.next_chain.branched_atom)
        self.assertEqual(rx.branched_atom, rx.next_chain.next_chain.branched_atom.ring_bonds[0].target)

    def test_parse_ring_bond(self):
        """Specific problematic ring bond cases (thus raising errors)"""

        wrong_smiles = [
            'C11',  # same atom
            'C12CCCCC12',  # same bond twice
            'C21CCCCC12',  # same bond twice
            'C12C2CCC1',  # "2" forms a direct pair
            'C-1CCCCC=1',  # not the same kind of bond
            'C12CCCCC1',  # unmatched "2"
        ]

        for s in wrong_smiles:
            with self.assertRaises(smiles_parser.ParserException, msg=s):
                smiles_parser.parse(s)

    def test_parse_chirality(self):

        correct_smiles = [
            ('[C@H](Br)(Cl)C', '@'),
            ('[C@TH1H](Br)(Cl)C', '@TH1'),  # = "@"
            ('[Xe@SP1H](F)(Br)Cl', '@SP1')
        ]

        for s, chirality in correct_smiles:
            p = smiles_parser.parse(s)
            self.assertEqual(p.chain.branched_atom.atom.chirality, chirality)

        self.assertEqual(
            smiles_parser.parse('C(Br)=[C@AL1]=C(O)C').chain.next_chain.branched_atom.atom.chirality, '@AL1')

        with self.assertRaises(smiles_parser.ParserException):
            smiles_parser.parse('[C@TH3](Br)(Cl)C')

        with self.assertRaises(smiles_parser.ParserException):
            smiles_parser.parse('[Xe@SP4](F)(Br)Cl')

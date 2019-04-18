from tests import OSmiPyTestCase

from osmipy import lexer, smiles_parser, smiles_ast


class ParserTestCase(OSmiPyTestCase):

    @staticmethod
    def parse(s):
        return smiles_parser.Parser(lexer.Lexer(s)).smiles()

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

        rx = smiles_parser_obj.smiles()

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

        rx = smiles_parser.Parser(lexer.Lexer('C1CC1')).smiles()

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

    def test_navigate_tree(self):
        """Test going down, up, left and right in the tree
        """

        smiles = ParserTestCase.parse('Cl\\C=C/F')

        # children
        children = list(smiles.children)

        self.assertIn(smiles.left, children)
        self.assertIn(smiles.bond, children)
        self.assertIn(smiles.right, children)

        self.assertNotIn(smiles.right.left, children)

        # contents
        self.assertEqual(children, smiles.contents)

        # descendants
        descendants = list(smiles.descendants)
        self.assertEqual(len(descendants), 14)

        self.assertEqual(len([i for i in descendants if type(i) is smiles_ast.Atom]), 4)  # Cl, C, C, F
        self.assertEqual(len([i for i in descendants if type(i) is smiles_ast.BranchedAtom]), 4)  # Same as above
        self.assertEqual(len([i for i in descendants if type(i) is smiles_ast.Bond]), 3)  # /, \, =
        self.assertEqual(len([i for i in descendants if type(i) is smiles_ast.Chain]), 3)

        # parent
        f_atom = descendants[-1]
        self.assertIs(type(f_atom), smiles_ast.Atom)
        self.assertIs(type(f_atom.parent), smiles_ast.BranchedAtom)

        # parents
        parents = list(f_atom.parents)
        self.assertEqual(len(parents), 5)  # BranchedAtom + 4 Chain
        self.assertEqual(len([i for i in parents if type(i) is smiles_ast.Chain]), 4)

        # next sibling(s)
        self.assertEqual(smiles.next_sibling, None)
        self.assertEqual(smiles.left.next_sibling, smiles.bond)
        self.assertEqual(smiles.bond.next_sibling, smiles.right)
        self.assertEqual(smiles.right.next_sibling, None)

        # previous sibling(s)
        self.assertEqual(smiles.previous_sibling, None)
        self.assertEqual(smiles.left.previous_sibling, None)
        self.assertEqual(smiles.bond.previous_sibling, smiles.left)
        self.assertEqual(smiles.right.previous_sibling, smiles.bond)

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
                smiles_parser.Parser(lexer.Lexer(s)).smiles()

    def test_hcount(self):
        """Test the specific case of hcount
        """

        to_test = [
            ('C', 4),
            ('c1ccccc1', 1),
            ('n1ccccc1', 0),  # /!\ should gives 1 (but aromatic detection is not yet implemented)
            ('N1C=CC=C1', 1),  # .. In the kekule form, that works
            ('N1C=CC=CC=1', 0),
            ('N(=O)O', 0),
            ('N(=O)(=O)O', 0),
            ('S(=O)=O', 0),
            ('SCCO', 1),
            ('C(=O)O', 1)
        ]

        for smi, hcount in to_test:
            node = ParserTestCase.parse(smi)
            self.assertEqual(hcount, node.left.implicit_hcount(), msg=smi)

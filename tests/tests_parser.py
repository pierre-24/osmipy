from tests import OSmiPyTestCase

from osmipy import lexer, smiles_parser, smiles_ast, smiles


class ParserTestCase(OSmiPyTestCase):

    @staticmethod
    def parse(s):
        """Parse a string to AST"""
        return smiles_parser.Parser(lexer.Lexer(s)).smiles()

    @staticmethod
    def interpret(ast):
        """Get a string representation of the AST"""
        return smiles.Interpreter(ast)()

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
        self.assertEqual(r.left.branches[0].chain.left.implicit_hcount, 0)
        self.assertEqual(r.bond.symbol, rx.bond.symbol)
        self.assertEqual(r.left.atom.symbol, rx.left.atom.symbol)  # sp² C
        self.assertEqual(r.left.implicit_hcount, 1)
        self.assertEqual(r.right.bond.symbol, rx.right.bond.symbol)
        self.assertEqual(r.right.left.atom.symbol, rx.right.left.atom.symbol)  # sp² C
        self.assertEqual(r.right.right.left.atom.symbol, rx.right.right.left.atom.symbol)  # -OH
        self.assertEqual(r.right.right.left.implicit_hcount, 1)

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
        self.assertEqual(rx.left.implicit_hcount, 2)  # sp³ C
        self.assertEqual(r.right.left.atom.symbol, rx.right.left.atom.symbol)
        self.assertEqual(rx.right.left.implicit_hcount, 2)  # sp³ C
        self.assertEqual(r.right.right.left.atom.symbol, rx.right.right.left.atom.symbol)
        self.assertEqual(rx.right.right.left.implicit_hcount, 2)  # sp³ C

        # ring bonds
        self.assertEqual(type(r.left.ring_bonds[0]), type(rx.left.ring_bonds[0]))
        self.assertEqual(type(r.right.right.left.ring_bonds[0]), type(rx.right.right.left.ring_bonds[0]))

        self.assertEqual(rx.left.ring_bonds[0].ring_id, rx.right.right.left.ring_bonds[0].ring_id)
        self.assertEqual(rx.left.ring_bonds[0].target, rx.right.right.left)
        self.assertEqual(rx.left, rx.right.right.left.ring_bonds[0].target)

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
            ('n1ccccc1', 0),
            ('N1C=CC=C1', 1),  # Note that this should be `[nH]1cccc1`, otherwise the implicit hcount is incorrect
            ('N1C=CC=CC=1', 0),
            ('N(=O)O', 0),
            ('N(=O)(=O)O', 0),
            ('S(=O)=O', 0),
            ('SCCO', 1),
            ('C(=O)O', 1)
        ]

        for smi, hcount in to_test:
            node = ParserTestCase.parse(smi)
            self.assertEqual(hcount, node.left.implicit_hcount, msg=smi)


class ASTTestCase(OSmiPyTestCase):

    def setUp(self):
        self.smiles = ParserTestCase.parse('Cl\\C=C/F')
        self.smiles2 = ParserTestCase.parse('C1(Cl)CC1')

    def test_navigate_tree(self):
        """Test going down, up, left and right in the tree
        """

        # children
        children = list(self.smiles.children)

        self.assertIn(self.smiles.left, children)
        self.assertIn(self.smiles.bond, children)
        self.assertIn(self.smiles.right, children)

        self.assertNotIn(self.smiles.right.left, children)

        # contents
        self.assertEqual(children, self.smiles.contents)

        # descendants
        descendants = list(self.smiles.descendants)
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
        self.assertEqual(self.smiles.next_sibling, None)
        self.assertEqual(self.smiles.left.next_sibling, self.smiles.bond)
        self.assertEqual(self.smiles.bond.next_sibling, self.smiles.right)
        self.assertEqual(self.smiles.right.next_sibling, None)

        # previous sibling(s)
        self.assertEqual(self.smiles.previous_sibling, None)
        self.assertEqual(self.smiles.left.previous_sibling, None)
        self.assertEqual(self.smiles.bond.previous_sibling, self.smiles.left)
        self.assertEqual(self.smiles.right.previous_sibling, self.smiles.bond)

    def test_remove(self):
        """Test node removal with `remove()`"""

        with self.assertRaises(ValueError):
            self.smiles.left.remove()  # cannot remove a `BranchedAtom`

        with self.assertRaises(ValueError):
            self.smiles.left.atom.remove()  # cannot remove an `Atom`

        with self.assertRaises(AttributeError):
            self.smiles.remove()

        # remove a bond in a chain
        self.assertIsNotNone(self.smiles.bond)
        self.smiles.bond.remove()
        self.assertIsNone(self.smiles.bond)

        # remove a chain
        self.assertIsNotNone(self.smiles.right.right.right)
        self.smiles.right.right.right.remove()
        self.assertIsNone(self.smiles.right.right.right)

        # remove a branch in a BranchedAtom
        self.assertNotEqual(self.smiles2.left.branches, [])
        self.smiles2.left.branches[0].remove()
        self.assertEqual(self.smiles2.left.branches, [])

        # remove a ring bond in a branched atom
        target = self.smiles2.left.ring_bonds[0].target
        self.assertNotEqual(self.smiles2.left.ring_bonds, [])
        self.smiles2.left.ring_bonds[0].remove()
        self.assertEqual(self.smiles2.left.ring_bonds, [])
        self.assertEqual(target.ring_bonds, [])  # both ends are deleted

        # indirectly delete a ring bond by removing the parent chain
        smiles3 = ParserTestCase.parse('C1CC1')
        self.assertNotEqual(smiles3.left.ring_bonds, [])
        smiles3.right.right.remove()
        self.assertEqual(smiles3.left.ring_bonds, [])

    def test_replace_with(self):
        """Test `replace_with()`
        """

        with self.assertRaises(AttributeError):  # cannot replace if no parent
            self.smiles.replace_with(ParserTestCase.parse('C'))

        with self.assertRaises(TypeError):  # cannot replace a node with another type of node
            self.smiles.right.replace_with(smiles_ast.Atom('C'))

        # replace atom
        self.assertEqual(self.smiles.left.atom.symbol, 'Cl')
        self.smiles.left.atom.replace_with(smiles_ast.Atom('H', isotope=2))  # change by a deuterium
        self.assertEqual(self.smiles.left.atom.symbol, 'H')  # atom changed
        self.assertEqual(self.smiles.left.atom.parent, self.smiles.left)  # parent is set !

        # replace bond on chain
        self.assertEqual(self.smiles.bond.symbol, '\\')
        self.smiles.bond.replace_with(smiles_ast.Bond('/'))
        self.assertEqual(self.smiles.bond.symbol, '/')

        # replace (right) chain
        self.smiles.right.right.right.replace_with(ParserTestCase.parse('OC'))
        self.assertEqual(self.smiles.right.right.right.left.atom.symbol, 'O')

        # replace branch
        self.assertEqual(len(list(self.smiles2.left.branches[0].descendants)), 3)  # Chain + BranchedAtom + Atom
        self.smiles2.left.branches[0].replace_with(smiles_ast.Branch(ParserTestCase.parse('COC')))
        self.assertNotEqual(len(list(self.smiles2.left.branches[0].descendants)), 3)

        # replace branched atom (which should also remove the other ring bond)
        target = self.smiles2.left.ring_bonds[0].target
        self.smiles2.left.replace_with(smiles_ast.BranchedAtom(atom=smiles_ast.Atom('O')))
        self.assertEqual(self.smiles2.left.ring_bonds, [])
        self.assertEqual(target.ring_bonds, [])  # both ends are deleted

    def test_insert(self):
        """Test `insert_before()` and `insert_after()` for branches and ring bonds
        """

        # test with branch:
        branch = smiles_ast.Branch(ParserTestCase.parse('CC'))
        self.assertNotEqual(self.smiles2.left.branches[0], branch)
        self.smiles2.left.branches[0].insert_before(branch)
        self.assertEqual(self.smiles2.left.branches[0], branch)

        branch.remove()  # remove itself from the AST
        self.smiles2.left.branches[0].insert_after(branch)
        self.assertEqual(self.smiles2.left.branches[1], branch)

        # test with ringbond
        ring_bond_beg = smiles_ast.RingBond(target=branch.chain.right.left, ring_id=2)
        ring_bond_end = smiles_ast.RingBond(target=self.smiles2.left, ring_id=2)

        branch.chain.right.left.append_ring_bond(ring_bond_end)

        self.assertNotEqual(self.smiles2.left.ring_bonds[0], ring_bond_beg)
        self.smiles2.left.ring_bonds[0].insert_before(ring_bond_beg)
        self.assertEqual(self.smiles2.left.ring_bonds[0], ring_bond_beg)

        ring_bond_beg.remove(signal_remove=False)  # do not remove the other ring bond !
        self.smiles2.left.ring_bonds[0].insert_after(ring_bond_beg)
        self.assertEqual(self.smiles2.left.ring_bonds[1], ring_bond_beg)

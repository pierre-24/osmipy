from tests import OSmiPyTestCase

from osmipy import smiles_parser, smiles_ast


class ASTTestCase(OSmiPyTestCase):

    def setUp(self):
        self.smiles = smiles_parser.parse('Cl\\C=C/F').chain
        self.smiles2 = smiles_parser.parse('C1(Cl)CC1').chain

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
            node = smiles_parser.parse(smi).chain
            self.assertEqual(hcount, node.left.implicit_hcount, msg=smi)

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
        self.assertEqual(len(parents), 6)  # BranchedAtom + 4 Chain + SMILES
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
        smiles3 = smiles_parser.parse('C1CC1').chain
        self.assertNotEqual(smiles3.left.ring_bonds, [])
        smiles3.right.right.remove()
        self.assertEqual(smiles3.left.ring_bonds, [])

    def test_replace_with(self):
        """Test `replace_with()`
        """

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
        self.smiles.right.right.right.replace_with(smiles_parser.parse('OC').chain)
        self.assertEqual(self.smiles.right.right.right.left.atom.symbol, 'O')

        # replace branch
        self.assertEqual(len(list(self.smiles2.left.branches[0].descendants)), 3)  # Chain + BranchedAtom + Atom
        self.smiles2.left.branches[0].replace_with(smiles_ast.Branch(smiles_parser.parse('COC').chain))
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
        branch = smiles_ast.Branch(smiles_parser.parse('CC').chain)
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

    def test_get_atom(self):
        s = smiles_parser.parse('CCO')
        self.assertEqual(s.get_atom(0), s.chain.left.atom)
        self.assertEqual(s.get_atom(1), s.chain.right.left.atom)
        self.assertEqual(s.get_atom(2), s.chain.right.right.left.atom)

    def test_add_fragment(self):
        """Test add fragment to a SMILES object"""

        a_str = 'c1ccccc1'
        b_str = 'N1CC2CCCC2CC1'

        a = smiles_parser.parse(a_str)
        next_atom_id = a._next_atom_id

        b = smiles_parser.parse(b_str)
        a += b

        self.assertEqual(str(a), a_str + '.' + b_str)
        self.assertEqual(a.get_atom(next_atom_id).symbol, 'N')  # ok, atoms_id updated !


class SMILESTestCase(OSmiPyTestCase):

    def test_smiles(self):
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
            'c1c2c3c4cc1.Br2.Cl3.Cl4',
            'c1ccccc1c1ccccc1'  # two times the same ring id
            'c1cc2cc3.c1cc2cc3',  # a very weird way of  generating an aromatic structure
        ]

        for smi in test_smiles:
            smile = smiles_parser.parse(smi)
            self.assertEqual(smi, str(smile))

from tests import OSmiPyTestCase

from osmipy import smiles


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
            smile = smiles.SMILES(smi)
            self.assertEqual(smi, repr(smile))

    def test_add_fragment(self):
        """Test add fragment to a SMILES object"""

        # test add_fragment
        a = smiles.SMILES('c1ccccc1')
        b = smiles.SMILES('N1CC2CCCC2CC1')
        smile = a.add_fragment(b)
        self.assertEqual(repr(smile), repr(a) + '.' + repr(b))

        self.assertEqual(smile.get_atom(a.next_atom_id).symbol, 'N')  # ok, atoms_id updated

    def test_get_atom(self):
        s = smiles.SMILES('CCO')
        self.assertEqual(s.get_atom(0), s.node.left.atom)
        self.assertEqual(s.get_atom(1), s.node.right.left.atom)
        self.assertEqual(s.get_atom(2), s.node.right.right.left.atom)

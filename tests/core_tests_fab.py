import unittest
from abpytools import Fab, ChainCollection
from operator import itemgetter
from urllib import request


abnum_url = 'http://www.bioinf.org.uk/abs/abnum'


# Helper functions
def check_connection(URL, timeout=5):
    try:
        request.urlopen(url=URL, timeout=timeout)
        return True
    except request.URLError:
        return False


class FabCore(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.light_chain_collection = ChainCollection(path='./tests/Data/chain_collection_1_light.json')
        cls.heavy_chain_collection = ChainCollection(path='./tests/Data/chain_collection_1_heavy.json')

        cls.light_chain_collection.load(verbose=False, show_progressbar=False)
        cls.heavy_chain_collection.load(verbose=False, show_progressbar=False)

        cls.heavy_chain = cls.heavy_chain_collection[0]
        cls.light_chain = cls.light_chain_collection[0]

    def test_Fab_input_exception_1(self):
        # cannot instantiate an 'empty' Fab object
        self.assertRaises(ValueError, Fab)

    def test_Fab_input_exception_2(self):
        self.assertRaises(ValueError, Fab, 0, 0)

    def test_Fab_input_exception_3(self):
        self.assertRaises(ValueError, Fab, [0], [0])

    def test_Fab_input_exception_4(self):
        # error is raised when one of the chains is not a chain object
        self.assertRaises(ValueError, Fab, self.heavy_chain, [0])

    def test_Fab_input_exception_5(self):
        self.assertRaises(ValueError, Fab, self.heavy_chain, self.heavy_chain)

    def test_Fab_input_exception_6(self):
        self.assertRaises(ValueError, Fab, self.light_chain, self.light_chain)

    def test_Fab_input_1(self):
        fab = Fab(heavy_chain=self.heavy_chain,  light_chain=self.light_chain, load=False)
        self.assertEqual(fab.name, 'ID1')

    def test_Fab_input_2(self):
        fab = Fab(heavy_chain=self.heavy_chain,  light_chain=self.light_chain, load=False, name='TEST')
        self.assertEqual(fab.name, 'TEST')

    def test_Fab_sequence(self):
        fab = Fab(heavy_chain=self.heavy_chain, light_chain=self.light_chain, load=False)
        self.assertEqual(fab.sequence, self.light_chain.sequence + self.heavy_chain.sequence)

    def test_Fab_input_internal_name_1(self):
        fab = Fab(heavy_chain=self.heavy_chain, light_chain=self.light_chain, load=False)
        self.assertEqual(fab._internal_light_name, '1QP1')

    def test_Fab_input_internal_name_2(self):
        fab = Fab(heavy_chain=self.heavy_chain, light_chain=self.light_chain, load=False)
        self.assertEqual(fab._internal_heavy_name, 'test')

    def test_Fab_molecular_weight_1(self):
        fab = Fab(heavy_chain=self.heavy_chain, light_chain=self.light_chain, load=False)
        self.assertAlmostEqual(fab.molecular_weight(monoisotopic=False), 31804.522078999988)

    def test_Fab_molecular_weight_2(self):
        fab = Fab(heavy_chain=self.heavy_chain, light_chain=self.light_chain, load=False)
        self.assertAlmostEqual(fab.molecular_weight(monoisotopic=True), 31823.98309999999)

    def test_Fab_ec_1(self):
        fab = Fab(heavy_chain=self.heavy_chain, light_chain=self.light_chain, load=False)
        self.assertAlmostEqual(fab.extinction_coefficient(reduced=False), 70080.0)

    def test_Fab_ec_2(self):
        fab = Fab(heavy_chain=self.heavy_chain, light_chain=self.light_chain, load=False)
        self.assertAlmostEqual(fab.extinction_coefficient(reduced=True), 69705.0)

    def test_Fab_ec_3(self):
        fab = Fab(heavy_chain=self.heavy_chain, light_chain=self.light_chain, load=False)
        self.assertAlmostEqual(fab.extinction_coefficient(reduced=False, normalise=True), 2.203460244613224)

    def test_Fab_Hmatrix_shape(self):
        fab = Fab(heavy_chain=self.heavy_chain, light_chain=self.light_chain, load=False)
        self.assertEqual(fab.hydrophobicity_matrix().shape, (296,))

    def test_Fab_Hmatrix_sum(self):
        fab = Fab(heavy_chain=self.heavy_chain, light_chain=self.light_chain, load=False)
        self.assertAlmostEqual(fab.hydrophobicity_matrix().sum(), -20.289999999999999)

    def test_Fab_charge_shape(self):
        fab = Fab(heavy_chain=self.heavy_chain, light_chain=self.light_chain, load=False)
        self.assertEqual(fab.charge().shape, (296,))

    def test_Fab_charge_sum(self):
        fab = Fab(heavy_chain=self.heavy_chain, light_chain=self.light_chain, load=False)
        self.assertAlmostEqual(fab.charge().sum(), -4.5403839845374332)

    def test_Fab_total_charge(self):
        fab = Fab(heavy_chain=self.heavy_chain, light_chain=self.light_chain, load=False)
        self.assertAlmostEqual(fab.total_charge(), -4.962297393494294)

    def test_Fab_numbering_table_1(self):
        fab = Fab(heavy_chain=self.heavy_chain, light_chain=self.light_chain, load=False)
        self.assertEqual(fab.numbering_table().shape, (1, 296))

    def test_Fab_numbering_table_2(self):
        fab = Fab(heavy_chain=self.heavy_chain, light_chain=self.light_chain, load=False)
        self.assertEqual(fab.numbering_table(chain='light').shape, (1, 138))

    def test_Fab_numbering_table_3(self):
        fab = Fab(heavy_chain=self.heavy_chain, light_chain=self.light_chain, load=False)
        self.assertEqual(fab.numbering_table(chain='heavy').shape, (1, 158))

    def test_Fab_numbering_table_4(self):
        fab = Fab(heavy_chain=self.heavy_chain, light_chain=self.light_chain, load=False)
        self.assertEqual(fab.numbering_table(region='CDR3', chain='heavy').loc['ID1'].values[10], 'L')

    def test_Fab_numbering_table_5(self):
        fab = Fab(heavy_chain=self.heavy_chain, light_chain=self.light_chain, load=False)
        self.assertEqual(fab.numbering_table(region='CDR3', chain='light').loc['ID1'].values[10], '-')

    def test_Fab_numbering_table_6(self):
        fab = Fab(heavy_chain=self.heavy_chain, light_chain=self.light_chain, load=False)
        self.assertEqual(fab.numbering_table(as_array=True, region='CDR3', chain='both')[25], 'L')

    def test_Fab_numbering_table_exception(self):
        fab = Fab(heavy_chain=self.heavy_chain, light_chain=self.light_chain, load=False)
        self.assertRaises(ValueError, fab.numbering_table, False, 'all', 'medium')

    def test_Fab_len(self):
        fab = Fab(heavy_chain=self.heavy_chain, light_chain=self.light_chain, load=False)
        self.assertEqual(len(fab), 291)

    def test_Fab_sclicing_exception(self):
        # can only slice fab with index 0 (light chain) or 1 (heavy chain)
        fab = Fab(heavy_chain=self.heavy_chain, light_chain=self.light_chain, load=False)
        self.assertRaises(ValueError, itemgetter([2]), fab)

    def test_Fab_sclicing_1(self):
        fab = Fab(heavy_chain=self.heavy_chain, light_chain=self.light_chain, load=False)
        self.assertEqual(fab[0].chain, 'light')

    def test_Fab_sclicing_2(self):
        fab = Fab(heavy_chain=self.heavy_chain, light_chain=self.light_chain, load=False)
        self.assertEqual(fab[1].chain, 'heavy')

    @unittest.skipUnless(check_connection(URL=abnum_url), 'No internet connection, skipping test.')
    def test_Fab_numbering_table(self):
        fab = Fab(heavy_chain=self.heavy_chain, light_chain=self.light_chain, load=False)
        self.assertEqual(fab.germline_identity['Average', 'Total'].values[0], 94.25)

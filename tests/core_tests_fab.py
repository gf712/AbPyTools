import unittest
from abpytools import Fab, ChainCollection


class FabCore(unittest.TestCase):

    def setUp(self):
        self.light_chain_collection = ChainCollection(path='./tests/chain_collection_1_light.json')
        self.heavy_chain_collection = ChainCollection(path='./tests/chain_collection_1_heavy.json')

        self.light_chain_collection.load()
        self.heavy_chain_collection.load()

        self.heavy_chain = self.heavy_chain_collection[0]
        self.light_chain = self.light_chain_collection[0]

    def test_Fab_input_exception_1(self):
        # cannot instantiate an 'empty' Fab object
        self.assertRaises(ValueError, Fab)

    def test_Fab_input_exception_2(self):
        self.assertRaises(ValueError, Fab, 0, 0)

    def test_Fab_input_exception_3(self):
        self.assertRaises(ValueError, Fab, [0], [0])

    def test_Fab_input_exception_4(self):
        self.assertRaises(ValueError, Fab, self.heavy_chain, self.heavy_chain)

    def test_Fab_input_exception_5(self):
        self.assertRaises(ValueError, Fab, self.light_chain, self.light_chain)

    def test_Fab_input(self):
        fab = Fab(heavy_chain=self.heavy_chain,  light_chain=self.light_chain, load=False)
        self.assertEqual(fab.name, 'ID1')

    def test_Fab_sequence(self):
        fab = Fab(heavy_chain=self.heavy_chain, light_chain=self.light_chain, load=False)
        self.assertEqual(fab.sequence, self.heavy_chain.sequence + self.light_chain.sequence)

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

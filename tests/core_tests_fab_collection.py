import unittest
from abpytools import FabCollection, ChainCollection, Fab
from urllib import request


abnum_url = 'http://www.bioinf.org.uk/abs/abnum'


# Helper functions
def check_connection(URL, timeout=5):
    try:
        request.urlopen(url=URL, timeout=timeout)
        return True
    except request.URLError:
        return False


class FabCollectionCore(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.light_chain_collection = ChainCollection(path='./tests/Data/chain_collection_light_2_sequences.json')
        cls.heavy_chain_collection = ChainCollection(path='./tests/Data/chain_collection_heavy_2_sequences.json')

        cls.light_chain_collection.load(verbose=False, show_progressbar=False)
        cls.heavy_chain_collection.load(verbose=False, show_progressbar=False)

        cls.heavy_chain_1 = cls.heavy_chain_collection[0]
        cls.light_chain_1 = cls.light_chain_collection[0]
        cls.heavy_chain_2 = cls.heavy_chain_collection[1]
        cls.light_chain_2 = cls.light_chain_collection[1]

        cls.fab_1 = Fab(heavy_chain=cls.heavy_chain_1, light_chain=cls.light_chain_1, name='Fab1')
        cls.fab_2 = Fab(heavy_chain=cls.heavy_chain_2, light_chain=cls.light_chain_2, name='Fab2')

    def test_FabCollection_input_exception_1(self):
        self.assertRaises(ValueError, FabCollection, None, 0, 0)

    def test_FabCollection_input_exception_2(self):
        self.assertRaises(ValueError, FabCollection, None, [0, 0], [0, 0])

    def test_FabCollection_input_exception_3(self):
        self.assertRaises(ValueError, FabCollection, None, self.heavy_chain_1, [0, 0])

    def test_FabCollection_input_exception_4(self):
        self.assertRaises(ValueError, FabCollection, None, [0, 0], self.light_chain_1)

    def test_FabCollection_input_exception_5(self):
        self.assertRaises(ValueError, FabCollection)

    def test_FabCollection_input_exception_6(self):
        # instantiation of FabCollection checks if there are as many heavy as light chains
        self.assertRaises(ValueError, FabCollection, None, [self.heavy_chain_1, self.heavy_chain_2],
                          [self.light_chain_1])

    def test_FabCollection_input_exception_7(self):
        # instantiation of FabCollection checks if there are as many names as fabs
        self.assertRaises(ValueError, FabCollection, [self.fab_1], None, None, ['foo', 'bar'])

    def test_FabCollection_input_exception_8(self):
        # instantiation of FabCollection checks if names is a list of strings
        self.assertRaises(ValueError, FabCollection, [self.fab_1], None, None, [1])

    def test_FabCollection_input_1(self):
        self.assertIsInstance(FabCollection(light_chains=self.light_chain_collection,
                                            heavy_chains=self.light_chain_collection),
                              FabCollection)

    def test_FabCollection_input_2(self):
        # check if input with a list with fab objects works
        self.assertIsInstance(FabCollection(fab=[self.fab_1, self.fab_2]),
                              FabCollection)

    def test_FabCollection_input_3(self):
        # check if input with a list of Chain objects
        self.assertIsInstance(FabCollection(heavy_chains=[self.heavy_chain_1, self.heavy_chain_2],
                                            light_chains=[self.light_chain_1, self.light_chain_2]),
                              FabCollection)

    def test_FabCollection_input4(self):
        light_chain_collection = ChainCollection(path='./tests/Data/chain_collection_light_2_sequences.json')
        heavy_chain_collection = ChainCollection(path='./tests/Data/chain_collection_heavy_2_sequences.json')
        fab_collection = FabCollection(light_chains=light_chain_collection,
                                       heavy_chains=heavy_chain_collection,
                                       names=['Fab1', 'Fab2'])
        self.assertEqual(fab_collection.numbering_table()['Light']['CDR3']['L89'].loc['Fab1'], 'Q')

    def test_FabCollection_MW(self):
        fab_collection = FabCollection(light_chains=self.light_chain_collection,
                                       heavy_chains=self.heavy_chain_collection)
        self.assertAlmostEqual(fab_collection.molecular_weight()[0], 25078.143331999992)

    def test_FabCollection_ec_1(self):
        fab_collection = FabCollection(light_chains=self.light_chain_collection,
                                       heavy_chains=self.heavy_chain_collection)
        self.assertEqual(fab_collection.extinction_coefficient()[0], 52870)

    def test_FabCollection_ec_2(self):
        fab_collection = FabCollection(light_chains=self.light_chain_collection,
                                       heavy_chains=self.heavy_chain_collection)
        self.assertAlmostEqual(fab_collection.extinction_coefficient(normalise=True)[0], 2.1082102969136987)

    def test_FabCollection_charge(self):
        fab_collection = FabCollection(light_chains=self.light_chain_collection,
                                       heavy_chains=self.heavy_chain_collection)
        self.assertAlmostEqual(fab_collection.charge().sum(), -2.0354649111988743)

    def test_FabCollection_total_charge(self):
        fab_collection = FabCollection(light_chains=self.light_chain_collection,
                                       heavy_chains=self.heavy_chain_collection)
        self.assertAlmostEqual(fab_collection.total_charge()[0], 3.3658691836339365)

    def test_FabCollection_name_1(self):
        fab_collection = FabCollection(light_chains=self.light_chain_collection,
                                       heavy_chains=self.heavy_chain_collection)
        self.assertEqual(fab_collection._heavy_chains.names[0], 'Seq1')

    def test_FabCollection_name_2(self):
        fab_collection = FabCollection(light_chains=self.light_chain_collection,
                                       heavy_chains=self.heavy_chain_collection)
        self.assertEqual(fab_collection._heavy_chains.names[1], 'Seq2')

    def test_FabCollection_name_3(self):
        fab_collection = FabCollection(light_chains=self.light_chain_collection,
                                       heavy_chains=self.heavy_chain_collection)
        self.assertEqual(fab_collection._light_chains.names[0], 'LightSeq1')

    def test_FabCollection_name_4(self):
        fab_collection = FabCollection(light_chains=self.light_chain_collection,
                                       heavy_chains=self.heavy_chain_collection)
        self.assertEqual(fab_collection._light_chains.names[1], 'LightSeq2')

    @unittest.skipUnless(check_connection(URL=abnum_url), 'No internet connection, skipping test.')
    def test_FabCollection_igblast_query_1(self):
        fab_collection = FabCollection(light_chains=self.light_chain_collection,
                                       heavy_chains=self.heavy_chain_collection)
        fab_collection.igblast_server_query()
        self.assertAlmostEqual(fab_collection.germline_identity['Heavy', 'Total'].iloc[0], 93.8)

    @unittest.skipUnless(check_connection(URL=abnum_url), 'No internet connection, skipping test.')
    def test_FabCollection_igblast_query_1(self):
        # check germline assignment
        fab_collection = FabCollection(light_chains=self.light_chain_collection,
                                       heavy_chains=self.heavy_chain_collection,
                                       names=['Fab1', 'Fab2'])
        fab_collection.igblast_server_query()
        self.assertEqual(fab_collection.germline['Heavy', 'Assignment']['Fab1'], 'IGHV1-8*01')

    @unittest.skipUnless(check_connection(URL=abnum_url), 'No internet connection, skipping test.')
    def test_FabCollection_igblast_query_1(self):
        # check germline assignment score
        fab_collection = FabCollection(light_chains=self.light_chain_collection,
                                       heavy_chains=self.heavy_chain_collection,
                                       names=['Fab1', 'Fab2'])
        fab_collection.igblast_server_query()
        self.assertAlmostEqual(fab_collection.germline['Heavy', 'Score'].loc['Fab1'], 1.82e-66)

    def test_FabCollection_HMatrix(self):
        fab_collection = FabCollection(light_chains=self.light_chain_collection,
                                       heavy_chains=self.heavy_chain_collection)
        self.assertEqual(fab_collection.hydrophobicity_matrix().shape, (2, 296))

    def test_FabCollection_names(self):
        fab_collection = FabCollection(light_chains=self.light_chain_collection,
                                       heavy_chains=self.heavy_chain_collection,
                                       names=['Fab1', 'Fab2'])
        self.assertEqual(fab_collection.names[0], 'Fab1')

    def test_FabCollection_n_ab(self):
        fab_collection = FabCollection(light_chains=self.light_chain_collection,
                                       heavy_chains=self.heavy_chain_collection)
        self.assertEqual(fab_collection.n_ab, 2)

    def test_FabCollection_len(self):
        fab_collection = FabCollection(light_chains=self.light_chain_collection,
                                       heavy_chains=self.heavy_chain_collection)
        self.assertEqual(len(fab_collection), 2)

    def test_FabCollection_germline_identity(self):
        fab_collection = FabCollection(light_chains=self.light_chain_collection,
                                       heavy_chains=self.heavy_chain_collection,
                                       names=['Fab1', 'Fab2'])
        self.assertAlmostEqual(fab_collection.germline_identity['Average']['CDR3'].loc['Fab1'], 78.55)

    def test_FabCollection_slicing_1(self):
        fab_collection = FabCollection(light_chains=self.light_chain_collection,
                                       heavy_chains=self.heavy_chain_collection,
                                       names=['Fab1', 'Fab2'])
        self.assertIsInstance(fab_collection[0], Fab)

    def test_FabCollection_slicing_2(self):
        fab_collection = FabCollection(light_chains=self.light_chain_collection,
                                       heavy_chains=self.heavy_chain_collection,
                                       names=['Fab1', 'Fab2'])
        self.assertEqual(fab_collection[0].name, 'Fab1')

    def test_FabCollection_slicing_3(self):
        fab_collection = FabCollection(light_chains=self.light_chain_collection,
                                       heavy_chains=self.heavy_chain_collection,
                                       names=['Fab1', 'Fab2'])
        self.assertEqual(fab_collection[[0, 1]].names, ['Fab1', 'Fab2'])

    def test_FabCollection_slicing_4(self):
        fab_collection = FabCollection(light_chains=self.light_chain_collection,
                                       heavy_chains=self.heavy_chain_collection,
                                       names=['Fab1', 'Fab2'])
        self.assertEqual(fab_collection[[0, 1]]._heavy_chains.names, ['Seq1', 'Seq2'])

    def test_FabCollection_slicing_5(self):
        fab_collection = FabCollection(light_chains=self.light_chain_collection,
                                       heavy_chains=self.heavy_chain_collection,
                                       names=['Fab1', 'Fab2'])
        self.assertEqual(fab_collection[[0, 1]]._light_chains.names, ['LightSeq1', 'LightSeq2'])

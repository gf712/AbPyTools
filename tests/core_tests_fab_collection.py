import unittest
from abpytools import FabCollection, ChainCollection
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


class FabCollectionCore(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.light_chain_collection = ChainCollection(path='./tests/Data/chain_collection_light_2_sequences.json')
        cls.heavy_chain_collection = ChainCollection(path='./tests/Data/chain_collection_heavy_2_sequences.json')

        cls.light_chain_collection.load(verbose=False, show_progressbar=False)
        cls.heavy_chain_collection.load(verbose=False, show_progressbar=False)

        cls.heavy_chain = cls.heavy_chain_collection[0]
        cls.light_chain = cls.light_chain_collection[0]

    def test_FabCollection_input_exception_1(self):
        self.assertRaises(ValueError, FabCollection, None, 0, 0)

    def test_FabCollection_input_exception_2(self):
        self.assertRaises(ValueError, FabCollection, None, [0, 0], [0, 0])

    def test_FabCollection_input_exception_3(self):
        self.assertRaises(ValueError, FabCollection, None, self.heavy_chain, [0, 0])

    def test_FabCollection_input_exception_4(self):
        self.assertRaises(ValueError, FabCollection, None, [0, 0], self.light_chain)

    def test_FabCollection_input_1(self):
        self.assertIsInstance(FabCollection(light_chains=self.light_chain_collection,
                                            heavy_chains=self.light_chain_collection),
                              FabCollection)

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
                                       heavy_chains=self.heavy_chain_collection)
        fab_collection.igblast_server_query()
        self.assertEqual(fab_collection.germline['Seq1 - LightSeq1']['Heavy'][0], 'IGHV1-8*01')

    @unittest.skipUnless(check_connection(URL=abnum_url), 'No internet connection, skipping test.')
    def test_FabCollection_igblast_query_1(self):
        # check germline assignment score
        fab_collection = FabCollection(light_chains=self.light_chain_collection,
                                       heavy_chains=self.heavy_chain_collection)
        fab_collection.igblast_server_query()
        self.assertEqual(fab_collection.germline['Seq1 - LightSeq1']['Heavy'][1], 1.82e-66)


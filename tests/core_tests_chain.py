import unittest
from abpytools import Chain, ChainCollection
from urllib import request


abnum_url = 'http://www.bioinf.org.uk/abs/abnum'


# Helper functions
def check_connection(URL, timeout=5):
    try:
        request.urlopen(url=URL, timeout=timeout)
        return True
    except request.URLError:
        return False


def read_sequence_from_single_sequence_fasta(path):
    with open(path, 'r') as f:
        data = f.readlines()[1]
    return data


class ChainCore(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.test_sequence = read_sequence_from_single_sequence_fasta('./tests/Data/chain_collection_fasta_test.fasta')
        cls.heavy_chain_collection_object = ChainCollection(path='./tests/Data/chain_collection_1_heavy.json')
        cls.heavy_chain_collection_object.load()
        cls.heavy_chain_object = cls.heavy_chain_collection_object[0]
        cls.light_chain_collection_object = ChainCollection(path='./tests/Data/chain_collection_1_light.json')
        cls.light_chain_collection_object.load()
        cls.light_chain_object = cls.light_chain_collection_object[0]

    def test_Chain_load(self):
        chain = Chain(sequence=self.test_sequence)
        chain.load()
        self.assertIsInstance(chain, Chain)

    def test_Chain_light_chain(self):
        self.assertEqual(self.light_chain_object.chain, 'light')

    def test_Chain_heavy_chain(self):
        self.assertEqual(self.heavy_chain_object.chain, 'heavy')

    def test_Chain_charge(self):
        self.assertAlmostEqual(self.heavy_chain_object.ab_charge().sum(), 1.7497642167513607)

    def test_Chain_numbering_table_np(self):
        self.assertEqual((self.heavy_chain_object.ab_numbering_table(as_array=True, region='CDR1') != '-').sum(), 7)

    def test_Chain_numbering_table_pd(self):
        self.assertEqual((self.heavy_chain_object.ab_numbering_table(as_array=False, 
                                                                     region='CDR1').values != '-').sum(), 7)

    def test_Chain_sequence_len(self):
        self.assertEqual(len(self.heavy_chain_object), 184)

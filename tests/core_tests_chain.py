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

    def setUp(self):
        self.test_sequence = read_sequence_from_single_sequence_fasta('./tests/chain_collection_fasta_test.fasta')
        self.chain_collection_object = ChainCollection(path='./tests/chain_collection_fasta_test.fasta')
        self.chain_collection_object.load()
        self.chain_object = self.chain_collection_object[0]

    def test_Chain_load(self):
        chain = Chain(sequence=self.test_sequence)
        chain.load()
        self.assertIsInstance(chain, Chain)

    def test_Chain_charge(self):
        self.assertAlmostEqual(self.chain_object.ab_charge().sum(), 1.7497642167513607)

    def test_Chain_numbering_table_np(self):
        self.assertEqual((self.chain_object.ab_numbering_table(as_array=True, region='CDR1') != '-').sum(), 7)

    def test_Chain_numbering_table_pd(self):
        self.assertEqual((self.chain_object.ab_numbering_table(as_array=False, region='CDR1').values != '-').sum(), 7)

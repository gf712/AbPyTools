import unittest
from abpytools import Chain


def read_sequence(path):
    with open(path, 'r') as f:
        data = f.readlines()[1]
    return data


class ChainCore(unittest.TestCase):

    def setUp(self):
        self.test_sequence = read_sequence('./tests/chain_collection_fasta_test.fasta')

    def test_Chain_load(self):
        chain = Chain(sequence=self.test_sequence)
        chain.load()
        self.assertIsInstance(chain, Chain)

    def test_Chain_charge(self):
        chain = Chain(sequence=self.test_sequence)
        chain.load()
        self.assertAlmostEqual(chain.ab_charge().sum(), 1.7497642167513607)

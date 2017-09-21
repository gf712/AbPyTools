import unittest
from abpytools.features import ChainDomains


class AbDomainTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.ab_file = './tests/Data/chain_collection_heavy_2_sequences.json'

    def test_ChainDomains_instantiation(self):
        chain_domain = ChainDomains(path=self.ab_file, load=True)
        self.assertIsInstance(chain_domain, ChainDomains)

    def test_ChainDomains_cdr_length(self):
        chain_domain = ChainDomains(path=self.ab_file, load=True)
        self.assertEqual(chain_domain.cdr_lengths()[0, 0], 7)

    def test_ChainDomains_cdr_sequence(self):
        chain_domain = ChainDomains(path=self.ab_file, load=True)
        self.assertEqual(chain_domain.cdr_sequences()['Seq1']['CDR3'], 'GLRYTRAGMIWG')

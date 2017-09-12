import unittest
from abpytools import ChainCollection, SequenceAlignment


class AntibodyCollectionCore(unittest.TestCase):

    def setUp(self):
        self.ab_collection_1 = ChainCollection(path='./tests/sequence_alignment_seq_1.fasta')
        self.ab_collection_2 = ChainCollection(path='./tests/sequence_alignment_seq_2.fasta')

        self.ab_collection_1.load()
        self.ab_collection_2.load()

    def test_needleman_wunsch_score(self):
        sa = SequenceAlignment(self.ab_collection_1[0], self.ab_collection_2, 'needleman_wunsch', 'BLOSUM62')
        self.assertEquals(sa.score[self.ab_collection_2.names[0]], 426)

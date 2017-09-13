import unittest
from abpytools import ChainCollection, SequenceAlignment


class SequenceAlignmentTests(unittest.TestCase):

    def setUp(self):
        self.ab_collection_1 = ChainCollection(path='./tests/sequence_alignment_seq_1.fasta')
        self.ab_collection_2 = ChainCollection(path='./tests/sequence_alignment_seq_2.fasta')

        self.ab_collection_1.load(show_progressbar=False, verbose=False)
        self.ab_collection_2.load(show_progressbar=False, verbose=False)

        self.seq2_aligned = self.read_sequence_from_file('./tests/BLOSUM62_aligned_sequence')

    def test_sequence_alignment_target(self):
        sa = SequenceAlignment(self.ab_collection_1[0], self.ab_collection_2, 'needleman_wunsch', 'BLOSUM62')
        self.assertEqual(sa.target_sequence, self.ab_collection_1[0].sequence)

    def test_needleman_wunsch_score_BLOSUM62(self):
        sa = SequenceAlignment(self.ab_collection_1[0], self.ab_collection_2, 'needleman_wunsch', 'BLOSUM62')
        sa.align_sequences()
        self.assertEqual(sa.score[self.ab_collection_2.names[0]], 426)

    def test_needleman_wunsch_score_BLOSUM80(self):
        sa = SequenceAlignment(self.ab_collection_1[0], self.ab_collection_2, 'needleman_wunsch', 'BLOSUM80')
        sa.align_sequences()
        self.assertEqual(sa.score[self.ab_collection_2.names[0]], 452)

    def test_needleman_wunsch_score_BLOSUM45(self):
        sa = SequenceAlignment(self.ab_collection_1[0], self.ab_collection_2, 'needleman_wunsch', 'BLOSUM45')
        sa.align_sequences()
        self.assertEqual(sa.score[self.ab_collection_2.names[0]], 513)

    def test_needleman_wunsch_aligned_sequences(self):
        sa = SequenceAlignment(self.ab_collection_1[0], self.ab_collection_2, 'needleman_wunsch', 'BLOSUM62')
        sa.align_sequences()
        self.assertEqual(sa.aligned_sequences[self.ab_collection_2.names[0]], self.seq2_aligned)

    def test_alignment_exception_1(self):
        # catch exception when subsitution matrix is not known
        self.assertRaises(ValueError, SequenceAlignment, self.ab_collection_1[0],
                          self.ab_collection_2, 'needleman_wunsch', 'foo')

    def test_alignment_exception_2(self):
        # catch exception when algorithm is not known
        sa = SequenceAlignment(self.ab_collection_1[0], self.ab_collection_2, 'foo', 'BLOSUM62')
        self.assertRaises(ValueError, sa.align_sequences)

    def test_alignment_indel_sign(self):
        # if indel is positive the user receives a warning and indel will be set to -indel
        sa = SequenceAlignment(self.ab_collection_1[0], self.ab_collection_2, 'needleman_wunsch', 'BLOSUM62')
        sa.align_sequences(indel=-10)
        sa_score_1 = sa.score
        sa = SequenceAlignment(self.ab_collection_1[0], self.ab_collection_2, 'needleman_wunsch', 'BLOSUM62')
        sa.align_sequences(indel=10)
        sa_score_2 = sa.score
        # since the second indel is positive it will be set to
        # its negative, which runs the algorithms
        # with the same params again and will get the same scores
        self.assertEqual(sa_score_1, sa_score_2)

    def read_sequence_from_file(self, file):
        with open(file, 'r') as f:
            data = f.readlines()[0]
        return data

import unittest
from abpytools import ChainCollection, SequenceAlignment


def read_sequence_from_file(file):
    with open(file, 'r') as f:
        data = f.readlines()[0]
    return data


class SequenceAlignmentTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.ab_collection_1 = ChainCollection.load_from_file(path='./tests/Data/sequence_alignment_seq_1.json',
                                                             show_progressbar=False, verbose=False)
        cls.ab_collection_2 = ChainCollection.load_from_file(path='./tests/Data/sequence_alignment_seq_2.json',
                                                             show_progressbar=False, verbose=False)

        cls.seq2_aligned = read_sequence_from_file('./tests/Data/BLOSUM62_aligned_sequence')

    def test_sequence_alignment_target(self):
        sa = SequenceAlignment(self.ab_collection_1[0], self.ab_collection_2, 'needleman_wunsch', 'BLOSUM62')
        self.assertEqual(sa.target_sequence, self.ab_collection_1[0].sequence)

    def test_needleman_wunsch_score_BLOSUMXX(self):
        test_cases = [
            ("BLOSUM45", 513),
            ("BLOSUM62", 426),
            ("BLOSUM80", 452)
        ]
        for x, output in test_cases:
            with self.subTest(name=x):
                sa = SequenceAlignment(self.ab_collection_1[0], self.ab_collection_2, 'needleman_wunsch', x)
                sa.align_sequences()
                self.assertEqual(sa.score[self.ab_collection_2.names[0]], output)

    def test_needleman_wunsch_aligned_sequences(self):
        sa = SequenceAlignment(self.ab_collection_1[0], self.ab_collection_2, 'needleman_wunsch', 'BLOSUM62')
        sa.align_sequences()
        self.assertEqual(sa.aligned_sequences[self.ab_collection_2.names[0]], self.seq2_aligned)

    def test_alignment_exception_1(self):
        # catch exception when substitution matrix is not known
        self.assertRaises(ValueError, SequenceAlignment, self.ab_collection_1[0],
                          self.ab_collection_2, 'needleman_wunsch', 'foo')

    def test_alignment_exception_2(self):
        # catch exception when algorithm is not known
        sa = SequenceAlignment(self.ab_collection_1[0], self.ab_collection_2, 'foo', 'BLOSUM62')
        self.assertRaises(ValueError, sa.align_sequences)

    def test_alignment_exception_3(self):
        # catch error when user tries to print alignment before calling .align_sequences()
        sa = SequenceAlignment(self.ab_collection_1[0], self.ab_collection_2, 'foo', 'BLOSUM62')
        self.assertRaises(ValueError, sa.print_aligned_sequences)

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

    def test_alignment_print_string(self):
        sa = SequenceAlignment(self.ab_collection_1[0], self.ab_collection_2, 'needleman_wunsch', 'BLOSUM62')
        sa.align_sequences()
        self.assertEqual(len(sa._aligned_sequences_string()), 3)

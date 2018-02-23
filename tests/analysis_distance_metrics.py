import unittest
from abpytools.analysis.distance_metrics import *


class DistanceMetricsTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.vector1 = [0, 1, 2, 3]
        cls.vector2 = [3, 2, 1, 0]

        cls.seq1 = 'Python'
        cls.seq2 = 'Peithen'

    def test_cosine_distance_1(self):
        self.assertAlmostEqual(cosine_distance(self.vector1, self.vector2),
                               1.2810446253588492)

    def test_cosine_distance_2(self):
        self.assertAlmostEqual(cosine_distance(self.vector1, self.vector1), 0)

    def test_cosine_similarity(self):
        self.assertAlmostEqual(cosine_similarity(self.vector1, self.vector2),
                               -0.2810446253588492)

    def test_hamming_distance_1(self):
        self.assertEqual(hamming_distance('foo', 'bar'), 3)

    def test_hamming_distance_2(self):
        self.assertRaises(ValueError, hamming_distance, 'foo', 'barr')

    def test_levenshtein_distance(self):
        self.assertEqual(levenshtein_distance(self.seq1, self.seq2), 3)

    def test_euclidean_distance(self):
        self.assertAlmostEqual(euclidean_distance(self.vector1, self.vector2),
                               4.4721359549995796)

    def test_manhattan_distance(self):
        self.assertEqual(manhattan_distance(self.vector1, self.vector2), 8)

import unittest
from abpytools import ChainCollection
from abpytools.features.composition import *


class SequenceCompositionTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.chain = ChainCollection.load_from_file(path='./tests/Data/chain_collection_1_light.json')
        cls.chain.load()
        cls.sequence = cls.chain.sequences[0]

    def test_chou(self):
        chou_pseudo_aa = chou_pseudo_aa_composition(self.sequence)
        self.assertAlmostEqual(chou_pseudo_aa[0][10], 1)
        self.assertAlmostEqual(chou_pseudo_aa[0][25], 504)
        self.assertAlmostEqual(chou_pseudo_aa[0][-1], 585.75)

    def test_aa_composition(self):
        composition = aa_composition(self.sequence)
        self.assertEqual(composition['F'], 3)

    def test_aa_frequency(self):
        composition = aa_frequency(self.sequence)
        self.assertAlmostEqual(composition['F'], 0.02803738317757, delta=10e-9)

    def test_distance_to_first(self):
        distance = distance_to_first(self.sequence)
        self.assertEqual(distance['W'], 34)

    def test_amino_acid_distribution(self):
        composition = aa_composition(self.sequence)
        distance = distance_to_first(self.sequence)
        distribution = aa_distribution(self.sequence, composition, distance)
        self.assertEqual(distribution['V'], 1290.5)

    def test_triad_method(self):
        triad = triad_method(self.sequence)
        self.assertEqual(len(triad[0]), 343)
        self.assertEqual(triad[0][255], 0.2)

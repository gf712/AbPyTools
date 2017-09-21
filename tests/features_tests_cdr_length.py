import unittest
from abpytools import CDRLength
import os


class CDRLengthTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.ab_file = './tests/Data/chain_collection_heavy_2_sequences.json'

    def test_CDRLength_instantiation(self):
        cdr_length = CDRLength(path=self.ab_file, load=True)
        self.assertIsInstance(cdr_length, CDRLength)

    def test_CDRLength_plot(self):
        # cannot test if plot works but can check if there is a file that is created
        cdr_length = CDRLength(path=self.ab_file, load=True)
        cdr_length.plot_cdr(save=True, plot_path='./tests', plot_name='TestCDR.png')
        self.assertTrue(os.path.isfile('./tests/TestCDR.png'), True)

    @classmethod
    def tearDownClass(cls):
        if os.path.isfile('./tests/TestCDR.png'):
            os.remove('./tests/TestCDR.png')

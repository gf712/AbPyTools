import unittest
from abpytools import Chain, ChainCollection
from abpytools.core.flags import *
from urllib import request
from parameterized import parameterized

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
        cls.heavy_chain_collection_object = ChainCollection.load_from_file(
            path='./tests/Data/chain_collection_1_heavy.json',
            show_progressbar=False, verbose=False)
        cls.heavy_chain_object = cls.heavy_chain_collection_object[0]
        cls.light_chain_collection_object = ChainCollection.load_from_file(
            path='./tests/Data/chain_collection_1_light.json',
            show_progressbar=False, verbose=False)
        cls.light_chain_object = cls.light_chain_collection_object[0]

    def test_Chain_load(self):
        chain = Chain(sequence=self.test_sequence)
        chain.load()
        self.assertIsInstance(chain, Chain)

    @parameterized.expand([
        ("light_chain", "light_chain_object", "light"),
        ("heavy_chain", "heavy_chain_object", "heavy")
    ])
    def test_Chain_chain(self, name, input, expected):
        self.assertEqual(getattr(self, input).chain, expected)

    @parameterized.expand([
        (f"{HYDROPHOBICITY_FLAGS.EW}_hydrophobicity_set", HYDROPHOBICITY_FLAGS.EW, -8.01),
        (f"{HYDROPHOBICITY_FLAGS.HH}_hydrophobicity_set", HYDROPHOBICITY_FLAGS.HH, 109.15),
        (f"{HYDROPHOBICITY_FLAGS.KD}_hydrophobicity_set", HYDROPHOBICITY_FLAGS.KD, -22.2),
        (f"{HYDROPHOBICITY_FLAGS.MF}_hydrophobicity_set", HYDROPHOBICITY_FLAGS.MF, 40.55),
        (f"{HYDROPHOBICITY_FLAGS.WW}_hydrophobicity_set", HYDROPHOBICITY_FLAGS.WW, -4.96)
    ])
    def test_Chain_hydrophobicity(self, name, input, expected):
        self.assertAlmostEqual(self.heavy_chain_object.ab_hydrophobicity_matrix(input).sum(), expected)

    def test_Chain_hydrophobicity_raise_error(self):
        self.assertRaises(ValueError, self.heavy_chain_object.ab_hydrophobicity_matrix, "TEST_I_DONT_EXIST")

    def test_Chain_hydrophobicity_not_loaded(self):
        test_sequence = read_sequence_from_single_sequence_fasta(
            './tests/Data/chain_collection_fasta_test.fasta')
        obj = Chain(sequence=test_sequence, name="test")
        self.assertAlmostEqual(obj.ab_hydrophobicity_matrix(HYDROPHOBICITY_FLAGS.EW).sum(), -8.01)

    @parameterized.expand([
        (f"{PI_FLAGS.EMBOSS}_pka_set", PI_FLAGS.EMBOSS, 1.9497227766014302),
        (f"{PI_FLAGS.DTASELECT}_pka_set", PI_FLAGS.DTASELECT, 1.9404670785960516),
        (f"{PI_FLAGS.SOLOMON}_pka_set", PI_FLAGS.SOLOMON, 1.7993251336784457),
        (f"{PI_FLAGS.SILLERO}_pka_set", PI_FLAGS.SILLERO, 2.0235864050746666),
        (f"{PI_FLAGS.RODWELL}_pka_set", PI_FLAGS.RODWELL, 1.8135202627308433),
        (f"{PI_FLAGS.WIKIPEDIA}_pka_set", PI_FLAGS.WIKIPEDIA, 1.7497642167513607),
        (f"{PI_FLAGS.LEHNINGER}_pka_set", PI_FLAGS.LEHNINGER, 1.8081589709221348),
        (f"{PI_FLAGS.GRIMSLEY}_pka_set", PI_FLAGS.GRIMSLEY, 0.5277757339608096)
    ])
    def test_Chain_charge(self, name, input, expected):
        self.assertAlmostEqual(self.heavy_chain_object.ab_charge(pka_database=input).sum(), expected)

    def test_Chain_charge_not_aligned(self):
        # this number is different from above since unaligned sequences can
        # include amino acids that are not found in the CDR and FR definitions
        self.assertAlmostEqual(self.heavy_chain_object.ab_charge(pka_database=PI_FLAGS.EMBOSS, align=False).sum(),
                               1.8015020954105214)

    @parameterized.expand([
        (f"{CHAIN_FLAGS.CDR1}", CHAIN_FLAGS.CDR1, 7),
        (f"{CHAIN_FLAGS.CDR2}", CHAIN_FLAGS.CDR2, 5),
        (f"{CHAIN_FLAGS.CDR3}", CHAIN_FLAGS.CDR3, 15),
        (f"{CHAIN_FLAGS.FR1}", CHAIN_FLAGS.FR1, 25),
        (f"{CHAIN_FLAGS.FR2}", CHAIN_FLAGS.FR2, 19),
        (f"{CHAIN_FLAGS.FR3}", CHAIN_FLAGS.FR3, 41),
        (f"{CHAIN_FLAGS.FR4}", CHAIN_FLAGS.FR4, 11),
    ])
    def test_Chain_numbering_table_np(self, name, input, expected):
        self.assertEqual((self.heavy_chain_object.ab_numbering_table(as_array=True, region=input) != '-').sum(),
                         expected)

    def test_Chain_numbering_table_pd(self):
        self.assertEqual((self.heavy_chain_object.ab_numbering_table(as_array=False,
                                                                     region='CDR1').values != '-').sum(), 7)

    @parameterized.expand([
        ("reduced_normalise", [True, True], 2.6041130777737154),
        ("notreduced_normalise", [False, True], 2.6165944479701),
        ("notreduced_notnormalise", [False, False], 52410),
        ("reduced_notnormalise", [True, False], 52160)
        ])
    def test_Chain_extinction_coefficient(self, name, input, expected):
        self.assertAlmostEqual(self.heavy_chain_object.ab_ec(reduced=input[0], normalise=input[1]), expected)

    def test_Chain_sequence_len(self):
        self.assertEqual(len(self.heavy_chain_object), 184)

    def test_Chain_set_name(self):
        heavy_chain_object = ChainCollection.load_from_file(
            path='./tests/Data/chain_collection_1_heavy.json',
            show_progressbar=False, verbose=False)
        heavy_chain_object[0].set_name("new_name")
        self.assertEqual(heavy_chain_object[0].name, "new_name")

    @unittest.skipUnless(check_connection(URL=abnum_url), 'No internet connection, skipping test.')
    def test_Chain_unnumbered_sequene(self):
        test_seq = Chain(sequence='TEST')
        test_seq.load()
        self.assertEqual(test_seq.status, 'Unnumbered')

    @parameterized.expand([
        ("chothia_numbering", NUMBERING_FLAGS.CHOTHIA, "H82A"),
        ("chothia_numbering", NUMBERING_FLAGS.CHOTHIA_EXT, "H80"),
        ("chothia_numbering", NUMBERING_FLAGS.KABAT, "H82A")
    ])
    @unittest.skipUnless(check_connection(URL=abnum_url), 'No internet connection, skipping test.')
    def test_numbering_scheme_alignment(self, name, input, expected):
        test = Chain(sequence=self.test_sequence, name="test", numbering_scheme=input)
        self.assertEqual(test.ab_numbering()[82], expected)

    @unittest.skipUnless(check_connection(URL=abnum_url), 'No internet connection, skipping test.')
    def test_numbering_scheme_alignment_unknow_scheme(self):
        test = Chain(sequence=self.test_sequence, name="test", numbering_scheme="MyNumberingScheme123")
        self.assertRaises(ValueError, test.ab_numbering)

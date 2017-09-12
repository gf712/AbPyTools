import unittest
from abpytools import ChainCollection, Chain
from urllib import request

url = 'http://www.bioinf.org.uk/abs/abnum'


def check_connection(URL, timeout=5):
    try:
        request.urlopen(url=URL, timeout=timeout)
        return True
    except request.URLError:
        return False


class AntibodyCollectionCore(unittest.TestCase):

    def setUp(self):
        self.antibody_collection_1_name = 'test'

    def test_ChainCollection_length_0(self):
        antibody_collection = ChainCollection()
        self.assertEqual(len(antibody_collection), 0)

    def test_ChainCollection_instance(self):
        antibody_collection = ChainCollection()
        self.assertIsNone(antibody_collection._path)

    def test_ChainCollection_chain(self):
        antibody_collection_1 = ChainCollection(path='./tests/chain_collection_1.json')
        antibody_collection_1.load(show_progressbar=False, verbose=False)
        self.assertEqual(antibody_collection_1.chain, 'heavy')

    def test_ChainCollection_Hmatrix_shape(self):
        antibody_collection_1 = ChainCollection(path='./tests/chain_collection_1.json')
        antibody_collection_1.load(show_progressbar=False, verbose=False)
        # if this fails it means that abysis has been updated
        self.assertEqual(antibody_collection_1.hydrophobicity_matrix().shape, (1, 158))

    def test_ChainCollection_sequence_length(self):
        antibody_collection_1 = ChainCollection(path='./tests/chain_collection_1.json')
        antibody_collection_1.load(show_progressbar=False, verbose=False)
        self.assertEqual(len(antibody_collection_1.sequences), 1)

    def test_ChainCollection_obj_length(self):
        antibody_collection_1 = ChainCollection(path='./tests/chain_collection_1.json')
        antibody_collection_1.load(show_progressbar=False, verbose=False)
        self.assertEqual(len(antibody_collection_1), 1)

    def test_ChainCollection_slicing_1_obj(self):
        antibody_collection_1 = ChainCollection(path='./tests/chain_collection_1.json')
        antibody_collection_1.load(show_progressbar=False, verbose=False)
        # if returning a single chain abpytools automatically creates a new Chain object
        self.assertIsInstance(antibody_collection_1[0], Chain)

    def test_ChainCollection_cdr_regions_part1(self):
        antibody_collection_1 = ChainCollection(path='./tests/chain_collection_1.json')
        antibody_collection_1.load(show_progressbar=False, verbose=False)
        self.assertCountEqual(antibody_collection_1.ab_region_index().keys(),
                              [self.antibody_collection_1_name])

    def test_ChainCollection_cdr_regions_part2(self):
        antibody_collection_1 = ChainCollection(path='./tests/chain_collection_1.json')
        antibody_collection_1.load(show_progressbar=False, verbose=False)
        self.assertCountEqual(antibody_collection_1.ab_region_index()[self.antibody_collection_1_name],
                              ['CDR', 'FR'])

    def test_ChainCollection_cdr_regions_part3_cdr(self):
        antibody_collection_1 = ChainCollection(path='./tests/chain_collection_1.json')
        antibody_collection_1.load(show_progressbar=False, verbose=False)
        self.assertCountEqual(antibody_collection_1.ab_region_index()[self.antibody_collection_1_name]['CDR'],
                              ['CDR1', 'CDR2', 'CDR3'])

    def test_ChainCollection_cdr_regions_part3_fr(self):
        antibody_collection_1 = ChainCollection(path='./tests/chain_collection_1.json')
        antibody_collection_1.load(show_progressbar=False, verbose=False)
        self.assertCountEqual(antibody_collection_1.ab_region_index()[self.antibody_collection_1_name]['FR'],
                              ['FR1', 'FR2', 'FR3', 'FR4'])

    def test_ChainCollection_charge(self):
        antibody_collection_1 = ChainCollection(path='./tests/chain_collection_1.json')
        antibody_collection_1.load(show_progressbar=False, verbose=False)
        self.assertAlmostEqual(antibody_collection_1.total_charge[self.antibody_collection_1_name], 1.3278508)

    def test_ChainCollection_igblast_parser_germline(self):
        antibody_collection_1 = ChainCollection(path='./tests/chain_collection_1.json')
        antibody_collection_1.load(show_progressbar=False, verbose=False)
        antibody_collection_1.imgt_local_query('tests/chain_collection_1_igblast.html')
        self.assertEqual(antibody_collection_1.germline_identity[self.antibody_collection_1_name][0], 'IGHV4-34*01')

    def test_ChainCollection_igblast_parser_germline_score(self):
        antibody_collection_1 = ChainCollection(path='./tests/chain_collection_1.json')
        antibody_collection_1.load(show_progressbar=False, verbose=False)
        antibody_collection_1.imgt_local_query('tests/chain_collection_1_igblast.html')
        self.assertEqual(antibody_collection_1.germline_identity[self.antibody_collection_1_name][1], 9.11e-69)

    def test_ChainCollection_(self):
        antibody_collection_1 = ChainCollection(path='./tests/chain_collection_1.json')
        antibody_collection_1.load(show_progressbar=False, verbose=False)
        self.assertIsInstance(antibody_collection_1.get_object('test'), Chain)

    @unittest.skipUnless(check_connection(URL=url), 'No internet connection, skipping test.')
    def test_Chain_abysis_parser(self):
        antibody = ChainCollection(path='./tests/chain_collection_fasta_test.fasta', numbering_scheme='chothia')
        antibody.load(verbose=False, show_progressbar=False)
        self.assertEqual(antibody.chain, 'heavy')

    @unittest.skipUnless(check_connection(URL=url), 'No internet connection, skipping test.')
    def test_Chain_abysis_parser_chothia(self):
        antibody = ChainCollection(path='./tests/chain_collection_fasta_test.fasta', numbering_scheme='chothia')
        antibody.load(verbose=False, show_progressbar=False)
        self.assertEqual(antibody.numbering_table(as_array=True)[0][-1], 'S')

    @unittest.skipUnless(check_connection(URL=url), 'No internet connection, skipping test.')
    def test_Chain_abysis_parser_kabat(self):
        antibody = ChainCollection(path='./tests/chain_collection_fasta_test.fasta', numbering_scheme='kabat')
        antibody.load(verbose=False, show_progressbar=False)
        self.assertEqual(antibody.numbering_table(as_array=True)[0][-1], '-')

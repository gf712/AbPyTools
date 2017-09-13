import unittest
from abpytools import ChainCollection, Chain
from urllib import request

abnum_url = 'http://www.bioinf.org.uk/abs/abnum'
imgt_url = 'https://www.ncbi.nlm.nih.gov/igblast/'


def check_connection(URL, timeout=5):
    try:
        request.urlopen(url=URL, timeout=timeout)
        return True
    except request.URLError:
        return False


class ChainCollectionCore(unittest.TestCase):

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

    def test_ChainCollection_n_ab(self):
        antibody_collection_1 = ChainCollection(path='./tests/chain_collection_1.json')
        antibody_collection_1.load(show_progressbar=False, verbose=False)
        self.assertEqual(antibody_collection_1.n_ab, 1)

    def test_ChainCollection_numbering_scheme(self):
        antibody_collection_1 = ChainCollection(path='./tests/chain_collection_1.json')
        antibody_collection_1.load(show_progressbar=False, verbose=False)
        self.assertEqual(antibody_collection_1.numbering_scheme, 'chothia')

    def test_ChainCollection_numbering_scheme_kabat(self):
        antibody_collection_1 = ChainCollection(path='./tests/chain_collection_1.json', numbering_scheme='kabat')
        antibody_collection_1.load(show_progressbar=False, verbose=False)
        self.assertEqual(antibody_collection_1.numbering_scheme, 'kabat')

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

    @unittest.skipUnless(check_connection(URL=imgt_url), 'No internet connection, skipping test.')
    def test_ChainCollection_imgt_server_query(self):
        antibody_collection_1 = ChainCollection(path='./tests/chain_collection_1.json')
        antibody_collection_1.load(show_progressbar=False, verbose=False)
        antibody_collection_1.imgt_server_query()
        self.assertEqual(antibody_collection_1.germline_identity[self.antibody_collection_1_name][0], 'IGHV4-34*01')

    @unittest.skipUnless(check_connection(URL=imgt_url), 'No internet connection, skipping test.')
    def test_ChainCollection_imgt_server_query(self):
        antibody_collection_1 = ChainCollection(path='./tests/chain_collection_1.json')
        antibody_collection_1.load(show_progressbar=False, verbose=False)
        antibody_collection_1.imgt_server_query()
        self.assertEqual(antibody_collection_1.germline_identity[self.antibody_collection_1_name][1], 9.11e-69)

    @unittest.skipUnless(check_connection(URL=imgt_url), 'No internet connection, skipping test.')
    def test_ChainCollection_imgt_server_query(self):
        antibody_collection_1 = ChainCollection(path='./tests/chain_collection_1.json')
        antibody_collection_1.load(show_progressbar=False, verbose=False)
        antibody_collection_1.imgt_server_query()
        self.assertEqual(antibody_collection_1.germline[self.antibody_collection_1_name]['Total'], 96.9)

    def test_ChainCollection_slicing(self):
        antibody_collection_1 = ChainCollection(path='./tests/chain_collection_1.json')
        antibody_collection_1.load(show_progressbar=False, verbose=False)
        self.assertIsInstance(antibody_collection_1.get_object('test'), Chain)

    @unittest.skipUnless(check_connection(URL=abnum_url), 'No internet connection, skipping test.')
    def test_Chain_abysis_parser(self):
        antibody = ChainCollection(path='./tests/chain_collection_fasta_test.fasta', numbering_scheme='chothia')
        antibody.load(verbose=False, show_progressbar=False)
        self.assertEqual(antibody.chain, 'heavy')

    @unittest.skipUnless(check_connection(URL=abnum_url), 'No internet connection, skipping test.')
    def test_Chain_abysis_parser_chothia(self):
        antibody = ChainCollection(path='./tests/chain_collection_fasta_test.fasta', numbering_scheme='chothia')
        antibody.load(verbose=False, show_progressbar=False)
        self.assertEqual(antibody.numbering_table(as_array=True)[0][-1], '-')

    @unittest.skipUnless(check_connection(URL=abnum_url), 'No internet connection, skipping test.')
    def test_Chain_abysis_parser_kabat(self):
        antibody = ChainCollection(path='./tests/chain_collection_fasta_test.fasta', numbering_scheme='kabat')
        antibody.load(verbose=False, show_progressbar=False)
        self.assertEqual(antibody.numbering_table(as_array=True)[0][-1], '-')

    def test_ChainCollection_numbering_tableDataFrame(self):
        antibody_collection_1 = ChainCollection(path='./tests/chain_collection_1.json')
        antibody_collection_1.load(show_progressbar=False, verbose=False)
        self.assertEqual(antibody_collection_1.numbering_table(as_array=False)['CDR1']['H32'].values[0], 'Y')

    def test_ChainCollection_numbering_table_shape_np(self):
        antibody_collection_1 = ChainCollection(path='./tests/chain_collection_1.json')
        antibody_collection_1.load(show_progressbar=False, verbose=False)
        self.assertEqual(antibody_collection_1.numbering_table(as_array=True).shape, (1,158))

    def test_ChainCollection_numbering_table_shape_pd(self):
        antibody_collection_1 = ChainCollection(path='./tests/chain_collection_1.json')
        antibody_collection_1.load(show_progressbar=False, verbose=False)
        self.assertEqual(antibody_collection_1.numbering_table(as_array=False).shape, (1,158))

    def test_ChainCollecction_numbering_table_region_pd(self):
        antibody_collection_1 = ChainCollection(path='./tests/chain_collection_1.json')
        antibody_collection_1.load(show_progressbar=False, verbose=False)
        self.assertEqual(
            antibody_collection_1.numbering_table(region='CDR1').loc[self.antibody_collection_1_name].values[-1], 'Y')

    def test_ChainCollection_numbering_table_region_np(self):
        antibody_collection_1 = ChainCollection(path='./tests/chain_collection_1.json')
        antibody_collection_1.load(show_progressbar=False, verbose=False)
        self.assertEqual(antibody_collection_1.numbering_table(as_array=True, region='CDR1')[0][-1], 'Y')

    def test_ChainCollection_numbering_table_fr_region(self):
        antibody_collection_1 = ChainCollection(path='./tests/chain_collection_1.json')
        antibody_collection_1.load(show_progressbar=False, verbose=False)
        self.assertEqual(antibody_collection_1.numbering_table(region='FR1').loc['test'].values[0], 'Q')

    def test_ChainCollection_molecular_weight(self):
        antibody_collection_1 = ChainCollection(path='./tests/chain_collection_1.json')
        antibody_collection_1.load(show_progressbar=False, verbose=False)
        self.assertAlmostEqual(antibody_collection_1.molecular_weights(monoisotopic=False)[0], 20029.85217699999)

    def test_ChainCollection_molecular_weight_monoisotopic(self):
        antibody_collection_1 = ChainCollection(path='./tests/chain_collection_1.json')
        antibody_collection_1.load(show_progressbar=False, verbose=False)
        self.assertAlmostEqual(antibody_collection_1.molecular_weights(monoisotopic=True)[0], 20042.1121)

    def test_ChainCollection_ec(self):
        antibody_collection_1 = ChainCollection(path='./tests/chain_collection_1.json')
        antibody_collection_1.load(show_progressbar=False, verbose=False)
        self.assertAlmostEqual(antibody_collection_1.extinction_coefficients(reduced=False)[0], 52410.0)

    def test_ChainCollection_ec_reduced(self):
        antibody_collection_1 = ChainCollection(path='./tests/chain_collection_1.json')
        antibody_collection_1.load(show_progressbar=False, verbose=False)
        self.assertAlmostEqual(antibody_collection_1.extinction_coefficients(reduced=True)[0], 52160.0)

    def test_ChainCollection_charge(self):
        antibody_collection_1 = ChainCollection(path='./tests/chain_collection_1.json')
        antibody_collection_1.load(show_progressbar=False, verbose=False)
        self.assertAlmostEqual(antibody_collection_1.charge.sum(), 1.7497642167513607)

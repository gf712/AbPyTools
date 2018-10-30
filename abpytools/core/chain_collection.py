from .chain import Chain
import numpy as np
import logging
from abpytools.utils import PythonConfig, Download
import json
import os
import pandas as pd
from .helper_functions import numbering_table_sequences, numbering_table_region, numbering_table_multiindex
from operator import itemgetter
from urllib import parse
from math import ceil
from .base import CollectionBase
from ..features.composition import *
from ..analysis.distance_metrics import *
from ..core.cache import Cache
from multiprocessing import Manager, Process
from inspect import signature
from .utils import (json_ChainCollection_formatter, pb2_ChainCollection_formatter, pb2_ChainCollection_parser,
                    fasta_ChainCollection_parser, json_ChainCollection_parser)
from .flags import *

# setting up debugging messages
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

ipython_config = PythonConfig()
if ipython_config.ipython_info == 'notebook':
    from tqdm import tqdm_notebook as tqdm  # pragma: no cover
else:
    from tqdm import tqdm

if BACKEND_FLAGS.HAS_PROTO:
    from abpytools.core.formats import ChainCollectionProto


class ChainCollection(CollectionBase):
    """
    Object containing Chain objects and to perform analysis on the ensemble.

    """

    def __init__(self, antibody_objects=None, load=True, **kwargs):
        """

        Args:
            antibody_objects:
            load:
            **kwargs:
        """

        if antibody_objects is None:
            self.antibody_objects = []
        else:
            if isinstance(antibody_objects, ChainCollection):
                antibody_objects = antibody_objects.antibody_objects
            elif not isinstance(antibody_objects, list):
                raise ValueError("Expected a list, instead got object of type {}".format(type(antibody_objects)))

            elif not all(isinstance(obj, Chain) for obj in antibody_objects) and len(antibody_objects) > 0:
                raise ValueError("Expected a list containing objects of type Chain")

            self.antibody_objects = antibody_objects

            if len(set(x.numbering_scheme for x in antibody_objects)) == 1:
                self._numbering_scheme = antibody_objects[0].numbering_scheme
            else:
                raise ValueError("ChainCollection only support Chain objects with the same numbering scheme.")

            if len(set(x.chain for x in antibody_objects)) == 1:
                self._chain = antibody_objects[0].chain
            elif len(set(x.chain for x in antibody_objects)) == 0:
                self._chain = ''
            else:
                raise ValueError("ChainCollection only support Chain objects with the same chain type.")

            if load:
                self.load(**kwargs)

    def load(self, show_progressbar=True, n_threads=4, verbose=True):
        self.antibody_objects, self._chain = load_from_antibody_object(
            antibody_objects=self.antibody_objects,
            show_progressbar=show_progressbar,
            n_threads=n_threads, verbose=verbose)

    @classmethod
    def load_from_fasta(cls, path, numbering_scheme=NUMBERING_FLAGS.CHOTHIA, n_threads=20,
                        verbose=True, show_progressbar=True):
        if not os.path.isfile(path):
            raise ValueError("File does not exist!")
        with open(path, 'r') as f:
            antibody_objects = fasta_ChainCollection_parser(f, numbering_scheme=numbering_scheme)

        chain_collection = cls(antibody_objects=antibody_objects, load=True,
                               n_threads=n_threads, verbose=verbose,
                               show_progressbar=show_progressbar)

        return chain_collection

    @classmethod
    def load_from_pb2(cls, path, n_threads=20, verbose=True, show_progressbar=True):
        with open(path, 'rb') as f:
            proto_parser = ChainCollectionProto()
            proto_parser.ParseFromString(f.read())

        antibody_objects = pb2_ChainCollection_parser(proto_parser)

        chain_collection = cls(antibody_objects=antibody_objects, load=True,
                               n_threads=n_threads, verbose=verbose,
                               show_progressbar=show_progressbar)

        return chain_collection

    @classmethod
    def load_from_json(cls, path, n_threads=20, verbose=True, show_progressbar=True):

        with open(path, 'r') as f:
            data = json.load(f)

        antibody_objects = json_ChainCollection_parser(data)

        chain_collection = cls(antibody_objects=antibody_objects, load=True,
                               n_threads=n_threads, verbose=verbose,
                               show_progressbar=show_progressbar)

        return chain_collection

    def save_to_json(self, path, update=True):
        with open(os.path.join(path + '.json'), 'w') as f:
            data = json_ChainCollection_formatter(self.antibody_objects)
            json.dump(data, f, indent=2)

    def save_to_pb2(self, path, update=True):
        proto_parser = ChainCollectionProto()
        try:
            with open(os.path.join(path + '.pb2'), 'rb') as f:
                proto_parser.ParseFromString(f.read())
        except IOError:
            # print("Creating new file")
            pass

        pb2_ChainCollection_formatter(self.antibody_objects, proto_parser)

        with open(os.path.join(path + '.pb2'), 'wb') as f:
            f.write(proto_parser.SerializeToString())

    def save_to_fasta(self, path, update=True):
        with open(os.path.join(path + '.fasta'), 'w') as f:
            f.writelines(make_fasta(self.names, self.sequences))

    def molecular_weights(self, monoisotopic=False):

        """

        :param monoisotopic: bool whether to use monoisotopic values
        :return: list
        """

        return [x.ab_molecular_weight(monoisotopic=monoisotopic) for x in self.antibody_objects]

    def extinction_coefficients(self, extinction_coefficient_database='Standard', reduced=False):

        """

        :param extinction_coefficient_database: string with the name of the database to use
        :param reduced: bool whether to consider the cysteines to be reduced
        :return: list
        """

        return [x.ab_ec(extinction_coefficient_database=extinction_coefficient_database,
                        reduced=reduced) for x in self.antibody_objects]

    def hydrophobicity_matrix(self):

        if self._chain == CHAIN_FLAGS.HEAVY_CHAIN:
            num_columns = 158
        else:
            num_columns = 138
        abs_hydrophobicity_matrix = np.zeros((len(self.antibody_objects), num_columns))

        for row in range(abs_hydrophobicity_matrix.shape[0]):
            abs_hydrophobicity_matrix[row] = self.antibody_objects[row].hydrophobicity_matrix

        return abs_hydrophobicity_matrix

    def get_object(self, name=''):

        """

        :param name: str
        :return:
        """

        if name in self.names:
            index = self.names.index(name)
            return self[index]
        else:
            raise ValueError('Could not find sequence with specified name')

    def ab_region_index(self):

        """
        method to determine index of amino acids in CDR regions
        :return: dictionary with names as keys and each value is a dictionary with keys CDR and FR
        'CDR' entry contains dictionaries with CDR1, CDR2 and CDR3 regions
        'FR' entry contains dictionaries with FR1, FR2, FR3 and FR4 regions
        """
        return {x.name: {'CDR': x.ab_regions()[0], 'FR': x.ab_regions()[1]} for x in self.antibody_objects}

    def numbering_table(self, as_array=False, region='all'):

        region = numbering_table_region(region)

        table = np.row_stack(
            [x.ab_numbering_table(as_array=True, region=region) for x in self.antibody_objects])

        if as_array:
            return table

        else:
            # return the data as a pandas.DataFrame -> it's slower but looks nicer and makes it easier to get
            # the data of interest
            whole_sequence_dict, whole_sequence = numbering_table_sequences(region, self._numbering_scheme, self._chain)

            multi_index = numbering_table_multiindex(region=region,
                                                     whole_sequence_dict=whole_sequence_dict)

            # create the DataFrame and assign the columns and index names
            data = pd.DataFrame(data=table)
            data.columns = multi_index
            data.index = self.names
            return data

    def igblast_server_query(self, chunk_size=50, show_progressbar=True, **kwargs):
        """

        :param show_progressbar:
        :param chunk_size:
        :param kwargs: keyword arguments to pass to igblast_options
        :return:
        """

        # check if query is larger than 50 sequences
        # if so split into several queries
        query_list = self._split_to_chunks(chunk_size=chunk_size)
        n_chunks = ceil(len(self) / chunk_size) - 1

        if show_progressbar:
            for query in tqdm(query_list, total=n_chunks):
                self._igblast_server_query(query, **kwargs)

        else:
            for query in query_list:
                self._igblast_server_query(query, **kwargs)

    def _igblast_server_query(self, query, **kwargs):

        # prepare raw data
        fasta_query = make_fasta(names=query.names, sequences=query.sequences)

        # get url with igblast options
        url = igblast_options(sequences=fasta_query, **kwargs)

        # send and download query
        q = Download(url, verbose=False)

        try:
            q.download()
        except ValueError:  # pragma: no cover
            raise ValueError("Check the internet connection.")  # pragma: no cover

        igblast_result = q.html

        self._parse_igblast_query(igblast_result, query.names)

    def igblast_local_query(self, file_path):

        # load in file
        with open(file_path, 'r') as f:
            igblast_result = f.readlines()

        self._parse_igblast_query(igblast_result, self.names)

    def append(self, antibody_obj):

        self.antibody_objects += antibody_obj

    def pop(self, index=-1):

        if index > len(self):
            raise ValueError("The given index is outside the range of the object.")

        element_to_pop = self[index]

        self._destroy(index=index)

        return element_to_pop

    def _destroy(self, index):

        del self.antibody_objects[index]

    # def filter(self):
    #
    #     # TODO: complete method
    #     pass
    #
    def set_numbering_scheme(self, numbering_scheme, realign=True):
        if realign:
            try:
                self._numbering_scheme = numbering_scheme
                self.antibody_objects, self._chain = load_from_antibody_object(self.antibody_objects)
            except:
                print("Could not realign sequences, nothing has been changed.")
        else:
            self._numbering_scheme = numbering_scheme

    @property
    def names(self):
        return [x.name for x in self.antibody_objects]

    @property
    def sequences(self):
        return [x.sequence for x in self.antibody_objects]

    @property
    def aligned_sequences(self):
        return [x.aligned_sequence for x in self.antibody_objects]

    @property
    def n_ab(self):
        return len(self.sequences)

    @property
    def chain(self):
        if self._chain == '':
            chains = set([x.chain for x in self.antibody_objects])
            if len(chains) == 1:
                self._chain = next(iter(chains))
                return self._chain
            else:
                raise ValueError('Different types of chains found in collection!')
        else:
            return self._chain

    @property
    def numbering_scheme(self):
        return self._numbering_scheme

    @property
    def charge(self):
        return np.array([x.ab_charge() for x in self.antibody_objects])

    @property
    def total_charge(self):
        return {x.name: x.ab_total_charge() for x in self.antibody_objects}

    @property
    def germline_identity(self):
        return {x.name: x.germline_identity for x in self.antibody_objects}

    @property
    def germline(self):
        return {x.name: x.germline for x in self.antibody_objects}

    def _string_summary_basic(self):
        return "abpytools.ChainCollection Chain type: {}, Number of sequences: {}".format(self._chain,
                                                                                          len(self.antibody_objects))

    def __repr__(self):
        return "<%s at 0x%02x>" % (self._string_summary_basic(), id(self))

    def __len__(self):
        return len(self.antibody_objects)

    def __getitem__(self, indices):
        if isinstance(indices, int):
            return self.antibody_objects[indices]
        else:
            return ChainCollection(antibody_objects=list(itemgetter(*indices)(self.antibody_objects)))

    def __add__(self, other):
        if isinstance(other, ChainCollection):
            if self.numbering_scheme != other.numbering_scheme:
                raise ValueError("Concatenation requires ChainCollection "
                                 "objects to use the same numbering scheme.")
            else:
                new_object_list = self.antibody_objects + other.antibody_objects

        elif isinstance(other, Chain):
            if self.numbering_scheme != other.numbering_scheme:
                raise ValueError("Concatenation requires Chain object to use "
                                 "the same numbering scheme as ChainCollection.")
            else:
                new_object_list = self.antibody_objects + [other]

        else:
            raise ValueError("Concatenation requires other to be of type "
                             "ChainCollection, got {} instead".format(type(other)))
        return ChainCollection(antibody_objects=new_object_list, load=False)

    def _split_to_chunks(self, chunk_size=50):
        """
        Helper function to split ChainCollection into size chunk_size and returns generator
        :param chunk_size: int, size of each chunk
        :return: generator to iterate of each chunk of size chunk_size
        """

        if self.n_ab > chunk_size:

            for x in range(0, self.n_ab, chunk_size):
                yield self[range(x, min(x + chunk_size, self.n_ab))]

        else:
            yield self

    def _parse_igblast_query(self, igblast_result, names):

        igblast_result_dict = load_igblast_query(igblast_result, names)

        # unpack results
        for name in names:
            obj_i = self.get_object(name=name)
            obj_i.germline = igblast_result_dict[name][1]
            obj_i.germline_identity = igblast_result_dict[name][0]

    def loading_status(self):
        return [x.status for x in self.antibody_objects]

    def composition(self, method='count'):
        """
        Amino acid composition of each sequence. Each resulting list is organised alphabetically (see composition.py)
        :param method:
        :return:
        """
        if method == 'count':
            return [order_seq(aa_composition(seq)) for seq in self.sequences]
        elif method == 'freq':
            return [order_seq(aa_frequency(seq)) for seq in self.sequences]
        elif method == 'chou':
            return chou_pseudo_aa_composition(self.sequences)
        elif method == 'triad':
            return triad_method(self.sequences)
        elif method == 'hydrophobicity':
            return self.hydrophobicity_matrix()
        elif method == 'volume':
            return side_chain_volume(self.sequences)
        else:
            raise ValueError("Unknown method")

    def distance_matrix(self, feature=None, metric='cosine_similarity', multiprocessing=False):

        """
        Returns the distance matrix using a given feature and distance metric
        :param feature: string with the name of the feature to use
        :param metric: string with the name of the metric to use
        :param multiprocessing: bool to turn multiprocessing on/off (True/False)
        :return: list of lists with distances between all sequences of len(data) with each list of len(data)
                 when i==j M_i,j = 0
        """

        if feature is None:
            transformed_data = self.sequences
        elif isinstance(feature, str):
            # in this case the features are calculated using a predefined featurisation method (see self.composition)
            transformed_data = self.composition(method=feature)

        elif isinstance(feature, list):
            # a user defined list with vectors
            if len(feature) != self.n_ab:
                raise ValueError("Expected a list of size {}, instead got {}.".format(self.n_ab, len(feature)))
            else:
                transformed_data = feature
        else:
            raise TypeError("Unexpected input for feature argument.")

        if metric == 'cosine_similarity':
            distances = self._run_distance_matrix(transformed_data, cosine_similarity, multiprocessing=multiprocessing)

        elif metric == 'cosine_distance':
            distances = self._run_distance_matrix(transformed_data, cosine_distance, multiprocessing=multiprocessing)

        elif metric == 'hamming_distance':
            # be careful hamming distance only works when all sequences have the same length
            distances = self._run_distance_matrix(transformed_data, hamming_distance, multiprocessing=multiprocessing)

        elif metric == 'levenshtein_distance':
            distances = self._run_distance_matrix(transformed_data, levenshtein_distance,
                                                  multiprocessing=multiprocessing)

        elif metric == 'euclidean_distance':
            distances = self._run_distance_matrix(transformed_data, euclidean_distance, multiprocessing=multiprocessing)

        elif metric == 'manhattan_distance':
            distances = self._run_distance_matrix(transformed_data, manhattan_distance, multiprocessing=multiprocessing)

        elif callable(metric):
            # user defined metric function
            user_function_signature = signature(metric)

            # number of params should be two, can take args with defaults though
            default_params = sum(['=' in x for x in user_function_signature.parameters])

            if len(user_function_signature.parameters) - default_params > 2:
                raise ValueError("Expected a function with two parameters")
            else:
                distances = self._run_distance_matrix(transformed_data, metric, multiprocessing=multiprocessing)

        else:
            raise ValueError("Unknown distance metric.")

        return distances

    def _run_distance_matrix(self, data, metric, multiprocessing=False):

        """
        Helper function to setup the calculation of each entry in the distance matrix
        :param data: list with all sequences
        :param metric: function that takes two string and calculates distance
        :param multiprocessing: bool to turn multiprocessing on/off (True/False)
        :return: list of lists with distances between all sequences of len(data) with each list of len(data)
                 when i==j M_i,j = 0
        """

        if multiprocessing:
            with Manager() as manager:
                cache = manager.dict()
                matrix = manager.dict()

                jobs = [Process(target=self._distance_matrix,
                                args=(data, i, metric, cache, matrix)) for i in range(len(data))]

                for j in jobs:
                    j.start()
                for j in jobs:
                    j.join()

                # order the data
                return [matrix[x] for x in range(len(data))]

        else:
            cache = Cache(max_cache_size=(len(data) * (len(data) - 1)) / 2)
            matrix = Cache(max_cache_size=len(data))

            for i in range(len(data)):
                cache.update(i, self._distance_matrix(data, i, metric, cache, matrix))

            return [matrix[x] for x in range(len(data))]

    @staticmethod
    def _distance_matrix(data, i, metric, cache, matrix):

        """
        Function to calculate distance from the ith sequence of the ith row to the remaining entries in the same row
        :param data: list with all sequences
        :param i: int that indicates the matrix row being processed
        :param metric: function that takes two string and calculates distance
        :param cache: either a Manager or Cache object to cache results
        :param matrix: either a Manager or Cache object to store results in a matrix
        :return: None
        """

        row = []
        seq_1 = data[i]
        for j, seq_2 in enumerate(data):

            if i == j:
                row.append(0)
                continue

            keys = ('{}-{}'.format(i, j), '{}-{}'.format(j, i))
            if keys[0] not in cache or keys[1] not in cache:
                cache['{}-{}'.format(i, j)] = metric(seq_1, seq_2)
            if keys[0] in cache:
                row.append(cache[keys[0]])
            elif keys[1] in cache:
                row.append(cache[keys[0]])
            else:
                raise ValueError("Bug in row {} and column {}".format(i, j))

        matrix[i] = row


def load_antibody_object(antibody_object):
    antibody_object.load()
    return antibody_object


def load_from_antibody_object(antibody_objects, show_progressbar=True, n_threads=20, verbose=True):
    """

    Args:
        antibody_objects (list):
        show_progressbar (bool):
        n_threads (int):
        verbose (bool):

    Returns:

    """

    if verbose:
        print("Loading in antibody objects")

    from queue import Queue
    import threading

    q = Queue()
    for i in range(n_threads):
        t = threading.Thread(target=worker, args=(q,))
        t.daemon = True
        t.start()

    if show_progressbar:
        for antibody_object in tqdm(antibody_objects):
            q.put(antibody_object)

    else:
        for antibody_object in antibody_objects:
            q.put(antibody_object)

    q.join()

    # if show_progressbar:
    #     aprun = parallelexecutor(use_bar='tqdm', n_jobs=n_jobs, timeout=timeout)
    # else:
    #     aprun = parallelexecutor(use_bar='None', n_jobs=n_jobs, timeout=timeout)
    #
    # # load in objects in parallel
    # antibody_objects = aprun(total=len(antibody_objects))(
    #     delayed(load_antibody_object)(obj) for obj in antibody_objects)

    status = [x.status for x in antibody_objects]
    failed = sum([1 if x == 'Not Loaded' or x == 'Failed' else 0 for x in status])

    # remove objects that did not load
    while 'Not Loaded' in status:
        i = status.index('Not Loaded')
        del antibody_objects[i], status[i]

    while 'Failed' in status:
        i = status.index('Failed')
        del antibody_objects[i], status[i]

    if verbose:
        print("Failed to load {} objects in list".format(failed))

    loaded_obj_chains = [x.chain for x in antibody_objects if x.status == 'Loaded']

    if len(set(loaded_obj_chains)) == 1:
        chain = loaded_obj_chains[0]
    else:
        raise ValueError("All sequences must be of the same chain type: Light or Heavy",
                         set([x.chain for x in loaded_obj_chains]))

    n_ab = len(loaded_obj_chains)

    if n_ab == 0:
        raise ValueError("Could not find any heavy or light chains in provided file or list of objects")

    return antibody_objects, chain


def load_igblast_query(igblast_result, names):
    """

    :param names:
    :param igblast_result:
    :return:
    """

    try:
        from bs4 import BeautifulSoup
    except ImportError:
        raise ImportError("Please install bs4 to parse the IGBLAST html file:"
                          "pip install beautifulsoup4")

    # instantiate BeautifulSoup object to make life easier with the html text!
    if isinstance(igblast_result, list):
        soup = BeautifulSoup(''.join(igblast_result), "lxml")
    else:
        soup = BeautifulSoup(igblast_result, "lxml")

    # get the results found in <div id="content"> and return the text as a string
    results = soup.find(attrs={'id': "content"}).text

    # get query names
    query = re.compile('Query: (.*)')
    query_ids = query.findall(results)

    # make sure that all the query names in query are in self.names
    if not set(names).issubset(set(query_ids)):
        raise ValueError('Make sure that you gave the same names in ChainCollection as you gave'
                         'in the query submitted to IGBLAST')

    # regular expression to get tabular data from each region
    all_queries = re.compile('(Query: .*?)\n\n\n\n', re.DOTALL)

    # parse the results with regex and get a list with each query data
    parsed_results = all_queries.findall(results)

    # regex to get the FR and CDR information for each string in parsed results
    region_finder = re.compile('^([CDR\d|FR\d|Total].*)', re.MULTILINE)

    result_dict = {}

    # iterate over each string in parsed result which contains the result for individual queries
    for query_result in parsed_results:

        # get query name and get the relevant object
        query_i = query.findall(query_result)[0]

        # check if the query being parsed is part of the object
        # (not all queries have to be part of the object, but the object names must be a subset of the queries)

        if query_i not in set(names):
            continue

        # list with CDR and FR info for query result
        region_info = region_finder.findall(query_result)

        # get the data from region info with dict comprehension
        germline_identity = {x.split()[0].split('-')[0]: float(x.split()[-1]) for x in region_info}

        # get the top germline assignment
        v_line_assignment = re.compile('V\s{}\t.*'.format(query_i))

        # the top germline assignment is at the top (index 0)
        germline_result = v_line_assignment.findall(results)[0].split()

        # store the germline assignment and the bit score in a tuple as the germline attribute of Chain
        germline = (germline_result[2], float(germline_result[-2]))

        result_dict[query_i] = (germline_identity, germline)

    return result_dict


def worker(q):
    while True:
        item = q.get()
        load_antibody_object(item)
        q.task_done()


def make_fasta(names, sequences):
    file_string = ''
    for name, sequence in zip(names, sequences):
        file_string += '>{}\n'.format(name)
        file_string += '{}\n'.format(sequence)

    return file_string


def igblast_options(sequences, domain='imgt',
                    germline_db_V='IG_DB/imgt.Homo_sapiens.V.f.orf.p',
                    germline_db_D='IG_DB/imgt.Homo_sapiens.D.f.orf',
                    germline_db_J='IG_DB / imgt.Homo_sapiens.J.f.orf',
                    num_alignments_V=1, num_alignments_D=1, num_alignments_J=1):
    values = {"queryseq": sequences,
              "germline_db_V": germline_db_V,
              "germline_db_D": germline_db_D,
              "germline_db_J": germline_db_J,
              "num_alignments_V": str(num_alignments_V),
              "num_alignments_D": str(num_alignments_D),
              "num_alignments_J": str(num_alignments_J),
              "outfmt": "7",
              "domain": domain,
              "program": "blastp"}

    url = "http://www.ncbi.nlm.nih.gov/igblast/igblast.cgi?"
    url += parse.urlencode(values)

    return url

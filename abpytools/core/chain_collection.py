from .chain import Chain
import numpy as np
import logging
from joblib import Parallel, delayed
from abpytools.utils import PythonConfig, Download
import json
import os
import pandas as pd
import re
from .helper_functions import numbering_table_sequences, numbering_table_region, numbering_table_multiindex
from operator import itemgetter
from urllib import parse
from math import ceil
from .base import CollectionBase
from ..features.composition import *
from ..analysis.distance_metrics import *
from ..core.cache import Cache
from multiprocessing import Manager, Process

# setting up debugging messages
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

ipython_config = PythonConfig()
if ipython_config.ipython_info == 'notebook':
    from tqdm import tqdm_notebook as tqdm
else:
    from tqdm import tqdm


class ChainCollection(CollectionBase):
    """
    Object containing Chain objects and to perform analysis on the ensemble.

    Methods:
        load: parses a FASTA file and creates an Chain object for each entry, or loads in Chain objects provided
        in a list
        name: returns a list with the names obtained from each Chain object's name attribute
        sequence: returns a list with the sequences of all antibodies in this ChainCollection object
        hydrophobicity_matrix: returns a numpy array with the hydrophobicity of shape (number of ab, 158/138)
        cdr_lengths: returns a numpy array with the CDR lengths of each antibody

    Attributes:
        chain: type of chain inferred form AbNum result
        n_ab: number of Chain objects

    """

    def __init__(self, antibody_objects=None, path=None, numbering_scheme='chothia'):

        """

        :type antibody_objects: list of Chain objects
              path:             string to a FASTA format file
        """

        if antibody_objects is None:
            antibody_objects = []

        if isinstance(antibody_objects, ChainCollection):
            antibody_objects = antibody_objects.antibody_objects
        elif not isinstance(antibody_objects, list):
            raise ValueError("Expected a list, instead got object of type {}".format(type(antibody_objects)))

        elif not all(isinstance(obj, Chain) for obj in antibody_objects) and len(antibody_objects) > 0:
            raise ValueError("Expected a list containing objects of type Chain")

        # check if path is a string
        elif not isinstance(path, str) and path is not None:
            raise ValueError("Path expected a str, instead got object of type {}".format(type(path)))

        elif isinstance(path, str) and path is not None:
            # check if path exists
            if not os.path.isfile(path):
                raise ValueError("The provided file path does not exist."
                                 "Expected a string containing the path to a FASTA or JSON format file")
            # check if it has the right extentions
            # TODO: additional checks to see if file has right format (not only right extension)
            elif not path.endswith('.json') and not path.endswith('.fasta'):
                raise ValueError("Expected the path to a FASTA or JSON format file")

        self.antibody_objects = antibody_objects

        if len(self.antibody_objects) > 0:
            if self.antibody_objects[0].chain != '':
                self._numbering_scheme = antibody_objects[0].numbering_scheme
                self._chain = antibody_objects[0].chain
            else:
                self._chain = ''
                self._numbering_scheme = numbering_scheme
        else:
            self._chain = ''
            self._numbering_scheme = numbering_scheme

        self._path = path

    def load(self, n_threads=20, verbose=True, show_progressbar=True, **kwargs):

        names = list()
        if self._path is not None:

            if self._path.endswith(".json"):
                with open(self._path, 'r') as f:
                    data = json.load(f)

                    ordered_names = data.pop('ordered_names')

                    for key_i in ordered_names:
                        antibody_dict_i = data[key_i]
                        antibody_i = Chain(name=key_i, sequence=antibody_dict_i['sequence'])
                        antibody_i.numbering = antibody_dict_i['numbering']
                        antibody_i._chain = antibody_dict_i['chain']
                        antibody_i.mw = antibody_dict_i['MW']
                        antibody_i.CDR = antibody_dict_i["CDR"]
                        antibody_i._numbering_scheme = antibody_dict_i["numbering_scheme"]
                        antibody_i.pI = antibody_dict_i["pI"]
                        antibody_i._loading_status = 'Loaded'
                        names.append(key_i)
                        self.antibody_objects.append(antibody_i)

            else:
                with open(self._path, 'r') as f:
                    names = list()
                    sequences = list()
                    for line in f:
                        if line.startswith(">"):
                            names.append(line.replace("\n", "")[1:])
                        # if line is empty skip line
                        elif line.isspace():
                            pass
                        else:
                            sequences.append(line.replace("\n", ""))
                    if len(names) != len(sequences):
                        raise ValueError("Error reading file: make sure it is FASTA format")

                    self.antibody_objects = list()

                    for name, sequence in zip(names, sequences):
                        self.antibody_objects.append(Chain(name=name, sequence=sequence))

        self.antibody_objects, self._chain = load_from_antibody_object(
            antibody_objects=self.antibody_objects,
            show_progressbar=show_progressbar,
            n_threads=n_threads, verbose=verbose, **kwargs)

    def save(self, file_format='FASTA', file_path='./', file_name='Ab_collection', information='all'):

        if file_format == 'FASTA':
            with open(os.path.join(file_path, file_name + '.fasta'), 'w') as f:
                f.writelines(make_fasta(self.names, self.sequences))

        if file_format == 'json':
            if information == 'all':

                with open(os.path.join(file_path, file_name + '.json'), 'w') as f:
                    # if antibody does not have name, generate name:
                    # ID_chain_idi, where chain is heavy/light, idi is i = [1,..,N]
                    idi = 1
                    data = dict()
                    ordered_names = []
                    for antibody in self.antibody_objects:

                        ordered_names.append(antibody.name)

                        antibody_dict = antibody.ab_format()
                        if len(antibody_dict['name']) > 0:
                            key_i = antibody_dict['name']
                        else:
                            key_i = "ID_{}_{}".format(antibody.chain, idi)
                            idi += 1
                        antibody_dict.pop("name")
                        data[key_i] = antibody_dict

                    data['ordered_names'] = ordered_names
                    json.dump(data, f, indent=2)

    def molecular_weights(self, monoisotopic=False):
        return [x.ab_molecular_weight(monoisotopic=monoisotopic) for x in self.antibody_objects]

    def extinction_coefficients(self, extinction_coefficient_database='Standard', reduced=False):
        return [x.ab_ec(extinction_coefficient_database=extinction_coefficient_database,
                        reduced=reduced) for x in self.antibody_objects]

    def hydrophobicity_matrix(self):

        if self._chain == 'heavy':
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
        except ValueError:
            raise ValueError("Check the internet connection.")

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

    @property
    def names(self):
        return [x.name for x in self.antibody_objects]

    @property
    def sequences(self):
        return [x.sequence for x in self.antibody_objects]

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
        return ChainCollection(antibody_objects=new_object_list)

    def _split_to_chunks(self, chunk_size=50):
        """
        Helper function to split ChainCollection into size chunk_size and returns generator
        :param chunks:
        :return:
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

        transformed_data = self.composition(method=feature)

        if metric == 'cosine_similarity':
            distances = self._run_distance_matrix(transformed_data, cosine_similarity, multiprocessing=multiprocessing)

        elif metric == 'cosine_distance':
            distances = self._run_distance_matrix(transformed_data, cosine_distance, multiprocessing=multiprocessing)

        elif metric == 'hamming_distance':
            # be careful hamming distance only works when all sequences have the same length
            distances = self._run_distance_matrix(self.sequences, hamming_distance, multiprocessing=multiprocessing)

        elif metric == 'levenshtein_distance':
            distances = self._run_distance_matrix(self.sequences, levenshtein_distance, multiprocessing=multiprocessing)

        elif metric == 'euclidean_distance':
            distances = self._run_distance_matrix(transformed_data, euclidean_distance, multiprocessing=multiprocessing)

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


# the following block of code was obtained from
# http://stackoverflow.com/questions/37804279/how-can-we-use-tqdm-in-a-parallel-execution-with-joblib
all_bar_funcs = {
    'tqdm': lambda args: lambda x: tqdm(x, **args),
    'None': lambda args: iter,
}


# def parallelexecutor(use_bar='tqdm', **joblib_args):
#     def aprun(bar=use_bar, **tq_args):
#         def tmp(op_iter):
#             if str(bar) in all_bar_funcs.keys():
#                 bar_func = all_bar_funcs[str(bar)](tq_args)
#             else:
#                 raise ValueError("Value %s not supported as bar type" % bar)
#             return Parallel(**joblib_args)(bar_func(op_iter))
#
#         return tmp
#
#     return aprun


def load_from_antibody_object(antibody_objects, show_progressbar=True, n_threads=20, verbose=True, timeout=5):
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
    loaded_obj_chains = [x.chain for x in antibody_objects if x.status != 'Not Loaded']

    failed = sum([1 if x == 'Not Loaded' else 0 for x in status])

    # remove objects that did not load
    while 'Not Loaded' in status:
        i = status.index('Not Loaded')
        del antibody_objects[i], status[i]

    if verbose:
        print("Failed to load {} objects in list".format(failed))

    if len(set(loaded_obj_chains)) == 1:
        chain = loaded_obj_chains[0]
    else:
        raise ValueError("All sequences must of the same chain type: Light or Heavy")

    n_ab = len(loaded_obj_chains)

    if n_ab == 0:
        raise IOError("Could not find any heavy or light chains in provided file or list of objects")

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

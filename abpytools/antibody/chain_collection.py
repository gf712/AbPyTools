from .chain import Chain
import numpy as np
import logging
from joblib import Parallel, delayed
from abpytools.utils import PythonConfig
import json
from os import path
import pandas as pd
import re
from .helper_functions import numbering_table_sequences, numbering_table_region, numbering_table_multiindex

# setting up debugging messages
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

ipython_config = PythonConfig()
ipython_config.get_ipython_info()
if ipython_config.backend == 'notebook':
    from tqdm import tqdm_notebook as tqdm
else:
    from tqdm import tqdm


class ChainCollection:
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

        if not isinstance(antibody_objects, list):
            raise IOError("Expected a list, instead got object of type {}".format(type(antibody_objects)))

        elif not all(isinstance(obj, Chain) for obj in antibody_objects) and len(antibody_objects) > 0:
            raise IOError("Expected a list containing objects of type Chain")

        elif not isinstance(path, str) and path is not None:
            raise IOError("Expected a string containing the path to a FASTA format file")

        self.antibody_objects = antibody_objects

        if len(self.antibody_objects) > 0:
            if self.antibody_objects[0].chain != '':
                self._numbering_scheme = antibody_objects[0].numbering_scheme
                self._chain = antibody_objects[0].chain
            else:
                self._chain = ''
                self._path = path
                self._numbering_scheme = numbering_scheme
        else:
            self._chain = ''
            self._path = path
            self._numbering_scheme = numbering_scheme

    def load(self, show_progressbar=True, n_jobs=-1):

        names = list()
        if self._path is not None:

            if self._path.endswith(".json"):
                with open(self._path, 'r') as f:
                    data = json.load(f)

                    for key_i in data.keys():
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
                        raise IOError("Error reading file: make sure it is FASTA format")

                    self.antibody_objects = list()

                    for name, sequence in zip(names, sequences):
                        self.antibody_objects.append(Chain(name=name, sequence=sequence))

        self.antibody_objects, self._chain = load_from_antibody_object(
            antibody_objects=self.antibody_objects,
            show_progressbar=show_progressbar,
            n_jobs=n_jobs)

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
        if name in self.names:
            index = self.names.index(name)
            return self.antibody_objects[index]
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

        table = np.array(
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

    def save(self, file_format='FASTA', file_path='./', file_name='Ab_FASTA.txt', information='all'):

        if file_format == 'FASTA':
            with open(path.join((file_path, file_name)), 'w') as f:
                for antibody in self.antibody_objects:
                    f.write('>{}\n'.format(antibody.name))
                    f.write('{}\n'.format(antibody.sequence))

        if file_format == 'json':
            if information == 'all':

                with open(path.join(file_path, file_name), 'w') as f:
                    # if antibody does not have name, generate name:
                    # ID_chain_idi, where chain is heavy/light, idi is i = [1,..,N]
                    idi = 1
                    data = dict()
                    for antibody in self.antibody_objects:

                        antibody_dict = antibody.ab_format()
                        if len(antibody_dict['name']) > 0:
                            key_i = antibody_dict['name']
                        else:
                            key_i = "ID_{}_{}".format(antibody.chain, idi)
                            idi += 1
                        antibody_dict.pop("name")
                        data[key_i] = antibody_dict
                    json.dump(data, f, indent=2)

    def append(self, antibody_obj):

        # TODO: complete method

        if isinstance(antibody_obj, Chain):
            self.antibody_objects.append(antibody_obj)
        elif isinstance(antibody_obj, ChainCollection):
            self.antibody_objects.append(antibody_obj.antibody_objects)

        self.antibody_objects = load_antibody_object(self.antibody_objects)

    def remove(self, antibody_obj=None, name=''):

        # TODO: complete method

        if isinstance(antibody_obj, Chain):
            name = antibody_obj.name
        elif isinstance(antibody_obj, ChainCollection):
            name = ChainCollection.names()
        elif antibody_obj is None and len(name) > 0:
            self.antibody_objects.pop([x for x in antibody_obj if x.name == name][0])

        self._update_obj()

    def filter(self):

        # TODO: complete method
        pass

    def _update_obj(self, index='all'):

        # TODO: write method
        if index == 'all':
            self.antibody_objects = load_from_antibody_object(self.antibody_objects, show_progressbar=True,
                                                              n_jobs=-1)

    def load_imgt_query(self, file_path):

        # TODO: Work on method to query IGBLAST directly using abpytools

        """
        Method to load in additional data from an IGBLAST file.
        Note for this to work properly you need to choose to get the data in tabular format and download the html
        Calling this method will give access to the germline identity of individual regions of the Fab

        :return: self
        """
        try:
            from bs4 import BeautifulSoup
        except ImportError:
            raise ImportError("Please install bs4 to parse the IGBLAST html file:"
                              "pip install beautifulsoup4")

        # load in file
        with open(file_path, 'r') as f:
            file = f.readlines()

        # instantiate BeautifulSoup object to make life easier with the html text!
        soup = BeautifulSoup(''.join(file), "lxml")

        # get the results found in <div id="content"> and return the text as a string
        results = soup.find(attrs={'id': "content"}).text

        # get query names
        query = re.compile('Query: (.*)')
        query_ids = query.findall(results)

        # make sure that all the query names in query are in self.names
        if not set(self.names).issubset(set(query_ids)):
            raise ValueError('Make sure that you gave the same names in ChainCollection as you gave'
                             'in the query submitted to IGBLAST')

        # regular expression to get tabular data from each region
        all_queries = re.compile('(Query: .*?)\n\n\n\n', re.DOTALL)

        # parse the results with regex and get a list with each query data
        parsed_results = all_queries.findall(results)

        # regex to get the FR and CDR information for each string in parsed results
        region_finder = re.compile('^([CDR\d|FR\d|Total].*)', re.MULTILINE)

        # iterate over each string in parsed result which contains the result for individual queries
        for query_result in parsed_results:

            # get query name and get the relevant object
            query_i = query.findall(query_result)[0]

            # check if the query being parsed is part of the object
            # (not all queries have to be part of the object, but the object names must be a subset of the queries)

            if query_i not in set(self.names):
                continue

            obj_i = self.get_object(name=query_i)

            # list with CDR and FR info for query result
            region_info = region_finder.findall(query_result)

            # get the data from region info with dict comprehension
            obj_i.germline_identity = {x.split()[0].split('-')[0]: float(x.split()[-1]) for x in region_info}

            # get the top germline assignment
            v_line_assignment = re.compile('V\s{}\t.*'.format(query_i))

            # the top germline assignment is at the top
            germline_result = v_line_assignment.findall(results)[0].split()

            # store the germline assignment and the bit score in a tuple as the germline attribute of Chain
            obj_i.germline = (germline_result[2], float(germline_result[-2]))

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


def load_antibody_object(antibody_object):
    antibody_object.load()
    return antibody_object


# the following block of code was obtained from
# http://stackoverflow.com/questions/37804279/how-can-we-use-tqdm-in-a-parallel-execution-with-joblib
all_bar_funcs = {
    'tqdm': lambda args: lambda x: tqdm(x, **args),
    'None': lambda args: iter,
}


def parallelexecutor(use_bar='tqdm', **joblib_args):
    def aprun(bar=use_bar, **tq_args):
        def tmp(op_iter):
            if str(bar) in all_bar_funcs.keys():
                bar_func = all_bar_funcs[str(bar)](tq_args)
            else:
                raise ValueError("Value %s not supported as bar type" % bar)
            return Parallel(**joblib_args)(bar_func(op_iter))

        return tmp

    return aprun


def load_from_antibody_object(antibody_objects, show_progressbar=True, n_jobs=-1):
    print("Loading in antibody objects")

    if show_progressbar:
        aprun = parallelexecutor(use_bar='tqdm', n_jobs=n_jobs)
    else:
        aprun = parallelexecutor(use_bar='None', n_jobs=n_jobs)

    # load in objects in parallel
    antibody_objects = aprun(total=len(antibody_objects))(
        delayed(load_antibody_object)(obj) for obj in antibody_objects)

    status = [x.status for x in antibody_objects]
    loaded_obj_chains = [x.chain for x in antibody_objects if x.status != 'Not Loaded']

    failed = sum([1 if x == 'Not Loaded' else 0 for x in status])

    # remove objects that did not load
    while 'Not Loaded' in status:
        i = status.index('Not Loaded')
        del antibody_objects[i], status[i]

    print("Failed to load {} objects in list".format(failed))

    if len(set(loaded_obj_chains)) == 1:
        chain = loaded_obj_chains[0]
    else:
        raise ValueError("All sequences must of the same chain type: Light or Heavy")

    n_ab = len(loaded_obj_chains)

    if n_ab == 0:
        raise IOError("Could not find any heavy or light chains in provided file or list of objects")

    return antibody_objects, chain

from .antibody import Antibody
import numpy as np
import logging
from joblib import Parallel, delayed
from abpytools.utils import PythonConfig
import json
from os import path
import pandas as pd
from ..utils import DataLoader
# setting up debugging messages
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)


ipython_config = PythonConfig()
ipython_config.get_ipython_info()
if ipython_config.backend == 'notebook':
    from tqdm import tqdm_notebook as tqdm
else:
    from tqdm import tqdm


class AntibodyCollection:
    """
    Object containing Antibody objects and to perform analysis on the ensemble.

    Methods:
        load: parses a FASTA file and creates an Antibody object for each entry, or loads in Antibody objects provided
        in a list
        name: returns a list with the names obtained from each Antibody object's name attribute
        sequence: returns a list with the sequences of all antibodies in this AntibodyCollection object
        hydrophobicity_matrix: returns a numpy array with the hydrophobicity of shape (number of ab, 158/138)
        cdr_lengths: returns a numpy array with the CDR lengths of each antibody

    Attributes:
        chain: type of chain inferred form AbNum result
        n_ab: number of Antibody objects

    """

    def __init__(self, antibody_objects=None, path=None):

        """

        :type antibody_objects: list of Antibody objects
              path:             string to a FASTA format file
        """

        if antibody_objects is None:
            antibody_objects = []

        if not isinstance(antibody_objects, list):
            raise IOError("Expected a list, instead got object of type {}".format(type(antibody_objects)))

        elif not all(isinstance(obj, Antibody) for obj in antibody_objects) and len(antibody_objects) > 0:
            raise IOError("Expected a list containing objects of type Antibody")

        elif not isinstance(path, str) and path is not None:
            raise IOError("Expected a string containing the path to a FASTA format file")

        self.antibody_objects = antibody_objects
        self.chain = ''
        self._path = path
        self.n_ab = 0

    def load(self, show_progressbar=True, n_jobs=-1, get_numbering=True):

        if self._path is not None:

            if self._path.endswith(".json"):
                with open(self._path, 'r') as f:
                    data = json.load(f)

                    for key_i in data.keys():
                        antibody_dict_i = data[key_i]
                        antibody_i = Antibody(name=key_i, sequence=antibody_dict_i['sequence'])
                        antibody_i.numbering = antibody_dict_i['numbering']
                        antibody_i.chain = antibody_dict_i['chain']
                        antibody_i.mw = antibody_dict_i['MW']
                        antibody_i.CDR = antibody_dict_i["CDR"]
                        antibody_i.numbering_scheme = antibody_dict_i["numbering_scheme"]
                        antibody_i.pI = antibody_dict_i["pI"]

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
                        self.antibody_objects.append(Antibody(name=name, sequence=sequence))

        if get_numbering:
            self.antibody_objects, self.chain, self.n_ab = load_from_antibody_object(
                antibody_objects=self.antibody_objects,
                show_progressbar=show_progressbar,
                n_jobs=n_jobs)

    def names(self):
        return [x.name for x in self.antibody_objects]

    def sequences(self):
        return [x.sequence for x in self.antibody_objects]

    def molecular_weights(self, monoisotopic=False):
        return [x.ab_molecular_weight(monoisotopic=monoisotopic) for x in self.antibody_objects]

    def extinction_coefficients(self, extinction_coefficient_database='Standard', reduced=False):
        return [x.ab_ec(extinction_coefficient_database=extinction_coefficient_database,
                        reduced=reduced) for x in self.antibody_objects]

    def hydrophobicity_matrix(self):

        if self.chain == 'heavy':
            num_columns = 158
        else:
            num_columns = 138
        abs_hydrophobicity_matrix = np.zeros((len(self.antibody_objects), num_columns))

        for row in range(abs_hydrophobicity_matrix.shape[0]):
            abs_hydrophobicity_matrix[row] = self.antibody_objects[row].hydrophobicity_matrix

        return abs_hydrophobicity_matrix

    def ab_region_index(self):

        """
        method to determine index of amino acids in CDR regions
        :return: dictionary with keys 'CDR' and 'Framework'
        'CDR' entry contains dictionaries with CDR1, CDR2 and CDR3 regions (keys) for each antibody in the same order as
        sequences
        'Framework' entry contains dictionaries with FR1, FR2, FR3 and FR4 regions
        """
        cdrs = [x.ab_regions()[0] for x in self.antibody_objects]
        frameworks = [x.ab_regions()[1] for x in self.antibody_objects]
        return {'CDRs': cdrs, 'Frameworks': frameworks}

    def numbering_table(self):

        idi = 1
        names = list()

        for antibody in self.antibody_objects:

            if len(antibody.name) > 0:
                names.append(antibody.name)
            else:
                names.append("ID_{}_{}".format(antibody.chain, idi))
                idi += 1

        data_loader = DataLoader(data_type='NumberingSchemes',
                                 data=[antibody.numbering_scheme, self.chain])
        whole_sequence_dict = data_loader.get_data()

        whole_sequence = whole_sequence_dict['withCDR']

        df = pd.DataFrame(columns=whole_sequence, index=names)

        for antibody, name in zip(self.antibody_objects, names):

            df.loc[name] = antibody.ab_numbering_table(name=name, only_array=True)

        return df

    def save(self, file_format='FASTA', file_path='./', file_name='Ab_FASTA.txt', information='all'):

        if file_format == 'FASTA':
            with open(path.join((file_path, file_name)), 'w') as f:
                for antibody in self.antibody_objects:
                    f.write('>{}\n'.format(antibody.name))
                    f.write('{}\n'.format(antibody.sequence))

        if file_format == 'json':
            if information == 'all':

                with open(path.join(file_path, file_name), 'w') as f:
                    # if antibody does not have name generate name:
                    # ID_chain_idi, where chain is heavy/light, idi is idi (1,..,N)
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

        if isinstance(antibody_obj, Antibody):
            self.antibody_objects.append(antibody_obj)
        elif isinstance(antibody_obj, AntibodyCollection):
            self.antibody_objects.append(antibody_obj.antibody_objects)

        self.antibody_objects = load_antibody_object(self.antibody_objects)

    def remove(self, antibody_obj=None, name=''):

        # TODO: complete method

        if isinstance(antibody_obj, Antibody):
            name = antibody_obj.name
        elif isinstance(antibody_obj, AntibodyCollection):
            name = AntibodyCollection.names()
        elif antibody_obj is None and len(name) > 0:
            self.antibody_objects.pop([x for x in antibody_obj if x.name == name][0])

        self._update_obj()

    def filter(self):
        pass

    def _update_obj(self, index='all'):

        # TODO: write method
        if index == 'all':
            self.antibody_objects == load_from_antibody_object(self.antibody_objects, show_progressbar=True,
                                                               n_jobs=-1)


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

    # updated from stackoverflow answer in
    # http://stackoverflow.com/questions/37804279/how-can-we-use-tqdm-in-a-parallel-execution-with-joblib
    if show_progressbar:
        aprun = parallelexecutor(use_bar='tqdm', n_jobs=n_jobs)
    else:
        aprun = parallelexecutor(use_bar='None', n_jobs=n_jobs)

    # load in objects in parallel
    antibody_objects = aprun(total=len(antibody_objects)) \
        (delayed(load_antibody_object)(obj) for obj in antibody_objects)

    chains = [x.chain for x in antibody_objects]
    chains_without_na = [x for x in chains if x != 'NA']

    skipped = len([x.chain for x in antibody_objects if x.chain == 'NA'])

    while 'NA' in chains:
        i = chains.index('NA')
        del antibody_objects[i], chains[i]

    print("Skipped {} objects in list".format(skipped))

    if len(set(chains_without_na)) == 1:
        chain = chains_without_na[0]
    else:
        raise ValueError("All sequences must of the same chain type: Light or Heavy")

    n_ab = len(chains_without_na)

    if n_ab == 0:
        raise IOError("Could not find any heavy or light chains in provided file or list of objects")

    return antibody_objects, chain, n_ab

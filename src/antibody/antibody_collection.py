from .antibody import Antibody
import numpy as np
import logging
from tqdm import tqdm
from joblib import Parallel, delayed

# setting up debugging messages
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)


class AntibodyCollection:

    def __init__(self, antibody_objects='', path=''):

        self._antibody_objects = antibody_objects
        self._chain = ''
        self.path = path

    def load_from_antibody_object(self, antibody_objects=None, show_progressbar=True, n_jobs=-1):

        print("Loading in antibody objects")

        if len(antibody_objects) > 0:

            with Parallel(n_jobs=n_jobs, verbose=5) as parallel:
                self._antibody_objects = parallel(delayed(load_antibody_object)(obj) for obj in antibody_objects)

            chains = [x.chain for x in self._antibody_objects]
            chains_without_na = [x for x in chains if x != 'NA']

            skipped = len([x.chain for x in self._antibody_objects if x.chain == 'NA'])

            while 'NA' in chains:
                i = chains.index('NA')
                del self._antibody_objects[i], chains[i]

            print("Skipped {} objects in list".format(skipped))

            if len(set(chains_without_na)) == 1:
                self._chain = chains_without_na[0]
            else:
                raise ValueError("All sequences must of the same chain type: Light or Heavy")

    def load_from_fasta(self, show_progressbar=True):

        with open(self.path, 'r') as f:
            names = []
            sequences = []
            for line in f:
                if line.startswith(">"):
                    names.append(line.replace("\n", "")[1:])
                else:
                    sequences.append(line.replace("\n", ""))
            if len(names) != len(sequences):
                raise IOError("Error reading file: make sure it is FASTA format")

        obj_list = []

        for name, sequence in zip(names, sequences):
            obj_list.append(Antibody(name=name, sequence=sequence))

        self.load_from_antibody_object(antibody_objects=obj_list, show_progressbar=show_progressbar)

    def names(self):
        return [x.name for x in self._antibody_objects]

    def sequences(self):
        return [x.sequence for x in self._antibody_objects]

    def hydrophobicity_matrix(self):

        if self._chain == 'heavy':
            num_columns = 158
        else:
            num_columns = 138
        abs_hydrophobicity_matrix = np.zeros((len(self._antibody_objects), num_columns))

        for row in range(abs_hydrophobicity_matrix.shape[0]):
            abs_hydrophobicity_matrix[row] = self._antibody_objects[row].hydrophobicity_matrix

        return abs_hydrophobicity_matrix


def load_antibody_object(antibody_object):
    antibody_object.load()
    return antibody_object

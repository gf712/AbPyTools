from .antibody import Antibody
import numpy as np
import logging
from tqdm import tqdm

# setting up debugging messages
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)


class AntibodyCollection:

    def __init__(self):

        self._antibody_objects = []
        self._chain = ''

    def load_from_antibody_object(self, antibody_objects=None, show_progressbar=True):

        print("Loading in antibody objects")

        if len(antibody_objects) > 0:
            chains = list()
            skipped = 0

            if show_progressbar:
                for antibody_object in tqdm(antibody_objects):
                    antibody_object.load()
                    self._antibody_objects.append(antibody_object)
                    if antibody_object.numbering == 'NA':
                        skipped += 1
                    else:
                        chains.append(antibody_object.chain)
            else:
                for antibody_object in antibody_objects:
                    antibody_object.load()
                    self._antibody_objects.append(antibody_object)
                    if antibody_object.numbering == 'NA':
                        skipped += 1
                    else:
                        chains.append(antibody_object.chain)

            print("Skipped {} objects in list".format(skipped))

            if len(set(chains)) == 1:
                self._chain = chains[0]
            else:
                raise ValueError("All sequences must of the same chain type: Light or Heavy")

    def load_from_fasta(self, path='', show_progressbar=True):

        with open(path, 'r') as f:
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

        if self._chain == 'Heavy':
            num_columns = 158
        else:
            num_columns = 138
        abs_hydrophobicity_matrix = np.zeros((len(self._antibody_objects), num_columns))

        for row in range(abs_hydrophobicity_matrix.shape[0]):
            abs_hydrophobicity_matrix[row] = self._antibody_objects[row].hydrophobicity_matrix

        return abs_hydrophobicity_matrix

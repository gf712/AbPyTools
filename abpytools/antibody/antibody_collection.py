from .antibody import Antibody
import numpy as np
import logging
from tqdm import tqdm
from joblib import Parallel, delayed

# setting up debugging messages
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)


class AntibodyCollection:
    def __init__(self, antibody_objects=None, path=None):

        """

        :type antibody_objects: list of Antibody objects
              path:             string to a FASTA format file
        """
        if antibody_objects is None:
            antibody_objects = []
        if not isinstance(antibody_objects, list):
            raise IOError("Expected a list, instead got object of type {}".format(type(antibody_objects)))
        if not all(isinstance(obj, Antibody) for obj in antibody_objects):
            raise IOError("Expected a list containing objects of type Antibody")
        if not isinstance(path, str):
            raise IOError("Expected a string containing the path to a FASTA format file")
        self._antibody_objects = antibody_objects
        self.chain = ''
        self._path = path
        self.n_ab = 0

    def load_from_antibody_object(self, show_progressbar=True, n_jobs=-1):

        print("Loading in antibody objects")

        # updated from stackoverflow answer in
        # http://stackoverflow.com/questions/37804279/how-can-we-use-tqdm-in-a-parallel-execution-with-joblib
        if show_progressbar:
            aprun = parallelexecutor(use_bar='tqdm', n_jobs=n_jobs)
        else:
            aprun = parallelexecutor(use_bar='None', n_jobs=n_jobs)
        self._antibody_objects = aprun(total=len(self._antibody_objects)) \
            (delayed(load_antibody_object)(obj) for obj in self._antibody_objects)

        chains = [x.chain for x in self._antibody_objects]
        chains_without_na = [x for x in chains if x != 'NA']

        skipped = len([x.chain for x in self._antibody_objects if x.chain == 'NA'])

        while 'NA' in chains:
            i = chains.index('NA')
            del self._antibody_objects[i], chains[i]

        print("Skipped {} objects in list".format(skipped))

        if len(set(chains_without_na)) == 1:
            self.chain = chains_without_na[0]
        else:
            raise ValueError("All sequences must of the same chain type: Light or Heavy")

        self.n_ab = len(chains_without_na)

        if self.n_ab == 0:
            raise IOError("Could not find any heavy or light chains in provided file or list of objects")

    def load_from_fasta(self, show_progressbar=True):

        with open(self._path, 'r') as f:
            names = []
            sequences = []
            for line in f:
                if line.startswith(">"):
                    names.append(line.replace("\n", "")[1:])
                else:
                    sequences.append(line.replace("\n", ""))
            if len(names) != len(sequences):
                raise IOError("Error reading file: make sure it is FASTA format")

        self._antibody_objects = list()

        for name, sequence in zip(names, sequences):
            self._antibody_objects.append(Antibody(name=name, sequence=sequence))

        self.load_from_antibody_object(show_progressbar=show_progressbar)

    def names(self):
        return [x.name for x in self._antibody_objects]

    def sequences(self):
        return [x.sequence for x in self._antibody_objects]

    def hydrophobicity_matrix(self):

        if self.chain == 'heavy':
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

from abpytools import ChainCollection
import numpy as np
from abpytools.core import Cache


class ChainDomains(ChainCollection):

    def __init__(self, antibody_objects=None, path=None, load=True, verbose=True, show_progressbar=True, n_threads=10):
        # expect a string which is a path to a FASTA file
        super().__init__(antibody_objects=antibody_objects, path=path)
        if load:
            self.load(verbose=verbose, show_progressbar=show_progressbar, n_threads=n_threads)

        self._cache = Cache(max_cache_size=5)

    def cdr_lengths(self):
        """
        method to obtain cdr_lengths
        :return: m by n matrix with CDR lengths, where m is the number of antibodies in ChainCollection and n is
        three, corresponding to the three CDRs.
        """

        if 'cdr_lengths' not in self._cache:

            cdr_length_matrix = np.zeros((self.n_ab, 3), dtype=np.int)
            cdr_sequences = self.cdr_sequences()

            for m, antibody in enumerate(self.antibody_objects):
                for n, cdr in enumerate(['CDR1', 'CDR2', 'CDR3']):
                    cdr_length_matrix[m, n] = len(cdr_sequences[antibody.name][cdr])

            self._cache.update(key='cdr_lengths', data=cdr_length_matrix)

        return self._cache['cdr_lengths']

    def cdr_sequences(self):
        """
        method that returns sequences of each cdr
        :return: list of dictionaries with keys 'CDR1', 'CDR2' and 'CDR3' containing a string with the respective amino
        acid sequence
        """

        if 'cdr_sequences' not in self._cache:

            cdr_sequences = dict()

            for antibody in self.antibody_objects:
                dict_i = dict()
                for cdr in ['CDR1', 'CDR2', 'CDR3']:
                    self.sequence_splitter_helper(antibody=antibody,
                                                  region=cdr,
                                                  index=0,
                                                  dict_i=dict_i)
                cdr_sequences[antibody.name] = dict_i

            self._cache.update(key='cdr_sequences', data=cdr_sequences)

        return self._cache['cdr_sequences']

    def framework_length(self):
        framework_length_matrix = np.zeros((self.n_ab, 4), dtype=np.int)
        fr_sequences = self.framework_sequences()
        for m, antibody in enumerate(self.antibody_objects):
            for n, framework in enumerate(['FR1', 'FR2', 'FR3', 'FR4']):
                framework_length_matrix[m, n] = len(fr_sequences[antibody.name][framework])

        return framework_length_matrix

    def framework_sequences(self):

        framework_sequences = dict()

        for antibody in self.antibody_objects:
            dict_i = dict()
            for framework in ['FR1', 'FR2', 'FR3', 'FR4']:
                self.sequence_splitter_helper(antibody=antibody,
                                              region=framework,
                                              index=1,
                                              dict_i=dict_i)
            framework_sequences[antibody.name] = dict_i

        return framework_sequences

    @staticmethod
    def sequence_splitter_helper(antibody, region, index, dict_i):
        seq_i = list()
        indices = antibody.ab_regions()[index][region]
        for i in indices:
            seq_i.append(antibody.sequence[i])
        dict_i[region] = ''.join(seq_i)

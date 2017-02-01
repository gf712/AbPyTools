from abpytools import AntibodyCollection
from abpytools.metrics import CDR
from collections import Counter
import numpy as np
from matplotlib import pyplot as plt
from abpytools.utils.data_loader import DataLoader

amino_acid_index = {"R": 0,
                    "N": 1,
                    "D": 2,
                    "E": 3,
                    "Q": 4,
                    "K": 5,
                    "S": 6,
                    "T": 7,
                    "C": 8,
                    "H": 9,
                    "M": 10,
                    "A": 11,
                    "V": 12,
                    "G": 13,
                    "I": 14,
                    "L": 15,
                    "F": 16,
                    "P": 17,
                    "W": 18,
                    "Y": 19}


class AminoAcidFreq:
    def __init__(self, antibodies, region='CDR3'):

        # expect a string which is a path to a FASTA file
        if isinstance(antibodies, str):
            self._antibodies = AntibodyCollection(path=antibodies)
            self._antibodies.load()
        # can also be a AntibodyCollection object
        elif isinstance(antibodies, AntibodyCollection):
            self._antibodies = antibodies
            # check if AntibodyCollection has been loaded (should have n_ab > 0)
            # TODO come up with a more elegant way to check if .load() method has been called
            if self._antibodies.n_ab == 0:
                self._antibodies.load()
        if region in ['all', 'CDRs', 'FRs', 'FR1', 'FR2', 'FR3', 'FR4', 'CDR1', 'CDR2', 'CDR3']:
            # get the sequence for the specified region
            self.region = region
            self._region_assignment = CDR(self._antibodies)
            if self.region.startswith('CDR'):
                self._sequences = [x[self.region] for x in self._region_assignment.cdr_sequences()]
                data_loader = DataLoader(data_type='CDR_positions', data=['chothia', self._antibodies.chain])
                self._numbering = data_loader.get_data()[self.region]
            elif self.region.startswith('FR'):
                self._sequences = [x[self.region] for x in self._region_assignment.framework_sequences()]
            elif self.region == 'all':
                self._sequences = [x.sequence for x in self._antibodies]
        else:
            raise ValueError()

        self._aa_freq = np.zeros((20, len(max(self._sequences, key=len))))
        self._aa_hyd_freq = np.zeros((3, len(max(self._sequences, key=len))))
        self._aa_chg_freq = np.zeros((3, len(max(self._sequences, key=len))))

    def _amino_acid_freq(self):

        for position in range(len(max(self._sequences, key=len))):

            # AVJHVILSVCIUV
            # AABIUBXIB
            # BAIDVSVCU

            position_sequence = [x[position] for x in self._sequences if len(x) > position]
            count_i = Counter(position_sequence)
            total_i = len(position_sequence)

            for amino_acid_i in count_i.keys():

                self._aa_freq[amino_acid_index[amino_acid_i], position] = count_i[amino_acid_i] / total_i

                # _aa_hyd_freq: row1 -> hydrophilic
                #               row2 -> moderate
                #               row3 -> hydrophobic
                if amino_acid_i in ['R', 'N', 'D', 'E', 'Q', 'K', 'S', 'T']:
                    self._aa_hyd_freq[0, position] += count_i[amino_acid_i]
                elif amino_acid_i in ['C', 'H', 'M']:
                    self._aa_hyd_freq[1, position] += count_i[amino_acid_i]
                else:
                    self._aa_hyd_freq[2, position] += count_i[amino_acid_i]

                # _aa_chg_freq: row1 -> negative
                #               row2 -> positive
                #               row3 -> neutral
                if amino_acid_i in ['D', 'E']:
                    self._aa_chg_freq[0, position] += count_i[amino_acid_i]
                elif amino_acid_i in ['R', 'K', 'H']:
                    self._aa_chg_freq[1, position] += count_i[amino_acid_i]
                else:
                    self._aa_chg_freq[2, position] += count_i[amino_acid_i]

            # normalize values
            self._aa_chg_freq[:, position] /= total_i
            self._aa_hyd_freq[:, position] /= total_i

    def plot(self, sort_by='name'):

        if sort_by not in ['name', 'hydropathy', 'charge']:
            raise ValueError("Argument for sort_by not valid. Valid arguments are name, hydrophobicity and charge")

        # calculate frequency matrices
        self._amino_acid_freq()

        plt.figure(figsize=(8, 8))
        for position in range(self._aa_freq.shape[1]):
            previous = 0
            # color = iter(cm.rainbow(np.linspace(0, 1, self._aa_freq.shape[1])))
            color = iter(plt.get_cmap('Vega20').colors)

            if sort_by == 'name':
                for i, amino_acid in enumerate(sorted(amino_acid_index.keys())):
                    c = next(color)
                    plt.bar(position, self._aa_freq[amino_acid_index[amino_acid], position], bottom=previous,
                            label=amino_acid, color=c)
                    previous += self._aa_freq[amino_acid_index[amino_acid], position]
            elif sort_by == 'hydropathy':
                for i, prop_i in enumerate(['Hydrophilic', 'Moderate', 'Hydrophobic']):
                    c = next(color)
                    plt.bar(position, self._aa_hyd_freq[i, position], bottom=previous, label=prop_i, color=c)
                    previous += self._aa_hyd_freq[i, position]
            else:
                for i, prop_i in enumerate(['Negative', 'Positive', 'Neutral']):
                    c = next(color)
                    plt.bar(position, self._aa_chg_freq[i, position], bottom=previous, label=prop_i, color=c)
                    previous += self._aa_chg_freq[i, position]

        plt.xticks(np.arange(self._aa_freq.shape[1]), self._numbering, rotation=60)
        plt.margins(0.02)
        plt.xlabel('Position', size=16)
        plt.ylabel('Frequency', size=16)
        if sort_by == 'name':
            plt.legend(sorted(amino_acid_index.keys()), loc='center left', bbox_to_anchor=(1, 0.5))
        elif sort_by == 'hydropathy':
            plt.legend(['Hydrophilic', 'Moderate', 'Hydrophobic'], loc='center left', bbox_to_anchor=(1, 0.5))
        else:
            plt.legend(['Negative', 'Positive', 'Neutral'], loc='center left', bbox_to_anchor=(1, 0.5))

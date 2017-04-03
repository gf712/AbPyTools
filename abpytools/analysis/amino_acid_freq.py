from abpytools import AntibodyCollection
from abpytools.metrics import CDR
from collections import Counter
import numpy as np
from matplotlib import pyplot as plt
from abpytools.utils.data_loader import DataLoader
from os import path
from abpytools.utils import PythonConfig

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

        regions = ['all', 'CDRs', 'FRs', 'FR1', 'FR2', 'FR3', 'FR4', 'CDR1', 'CDR2', 'CDR3']
        if region in regions:
            # get the sequence for the specified region
            self.region = region
            self._region_assignment = CDR(self._antibodies)

            if self.region.startswith('CDR'):
                self._sequences = [x[self.region] for x in self._region_assignment.cdr_sequences()]
                data_loader = DataLoader(data_type='CDR_positions', data=['chothia', self._antibodies.chain])
                self._numbering = data_loader.get_data()[self.region]

            elif self.region.startswith('FR'):
                self._sequences = [x[self.region] for x in self._region_assignment.framework_sequences()]

            # TODO: implement 'all'
            elif self.region == 'all':
                self._sequences = [x.sequence for x in self._antibodies]
        else:
            raise ValueError('Parameter region must be either: {}. Not {}'.format(' ,'.join(regions), region))

        self._sequence_count = len(max(self._sequences, key=len))

        self._aa_freq = np.zeros((20, self._sequence_count))
        self._aa_hyd_freq = np.zeros((3, self._sequence_count))
        self._aa_chg_freq = np.zeros((3, self._sequence_count))

        self._aa_count = np.zeros((20, self._sequence_count))
        self._aa_hyd_count = np.zeros((3, self._sequence_count))
        self._aa_chg_count = np.zeros((3, self._sequence_count))

    def _amino_acid_freq(self, normalize):

        # if the sum of self._aa_count is zero then the count has not been performed at this point
        if self._aa_count.sum() == 0:

            for position in range(len(max(self._sequences, key=len))):

                position_sequence = [x[position] for x in self._sequences if len(x) > position]
                count_i = Counter(position_sequence)
                total_i = len(position_sequence)

                for amino_acid_i in count_i.keys():

                    self._aa_count[amino_acid_index[amino_acid_i], position] = count_i[amino_acid_i]

                    # _aa_hyd_freq: row1 -> hydrophilic
                    #               row2 -> moderate
                    #               row3 -> hydrophobic
                    if amino_acid_i in ['R', 'N', 'D', 'E', 'Q', 'K', 'S', 'T']:
                        self._aa_hyd_count[0, position] += count_i[amino_acid_i]
                    elif amino_acid_i in ['C', 'H', 'M']:
                        self._aa_hyd_count[1, position] += count_i[amino_acid_i]
                    else:
                        self._aa_hyd_count[2, position] += count_i[amino_acid_i]

                    # _aa_chg_freq: row1 -> negative
                    #               row2 -> positive
                    #               row3 -> neutral
                    if amino_acid_i in ['D', 'E']:
                        self._aa_chg_count[0, position] += count_i[amino_acid_i]
                    elif amino_acid_i in ['R', 'K', 'H']:
                        self._aa_chg_count[1, position] += count_i[amino_acid_i]
                    else:
                        self._aa_chg_count[2, position] += count_i[amino_acid_i]

                # normalize values
                # doing it even when it is not required comes at a small computational cost
                # it would take longer if the user had to recalculate everything to have a count plot and then a
                # frequency plot
                self._aa_freq[:, position] = self._aa_count[:, position] / total_i
                self._aa_chg_freq[:, position] = self._aa_chg_count[:, position] / total_i
                self._aa_hyd_freq[:, position] = self._aa_hyd_count[:, position] / total_i

        if normalize:
            return self._aa_freq, self._aa_chg_freq, self._aa_hyd_freq
        else:
            return self._aa_count, self._aa_chg_count, self._aa_hyd_count

    def plot(self, sort_by='name', normalize=True, display_count=True, plot_path='./',
             plot_name='AminoAcidFrequency.png', notebook_plot=True):

        if sort_by not in ['name', 'hydropathy', 'charge']:
            raise ValueError("Argument for sort_by not valid. Valid arguments are name, hydrophobicity and charge")

        # get count/ freq matrices
        # to avoid writing more code than necessary the count and freq are stored in the same variable
        # since they will always be plotted independently
        aa, chg, hyd = self._amino_acid_freq(normalize=normalize)

        fig = plt.figure(1, figsize=(8, 8))
        ax = fig.add_subplot(111)

        for position in range(self._aa_freq.shape[1]):
            previous = 0
            # 20 distinct colors
            colors = ["#023fa5", "#7d87b9", "#bec1d4", "#d6bcc0", "#bb7784", "#8e063b", "#4a6fe3", "#8595e1",
                      "#b5bbe3", "#e6afb9", "#e07b91", "#d33f6a", "#11c638", "#8dd593", "#c6dec7", "#ead3c6",
                      "#f0b98d", "#ef9708", "#0fcfc0", "#9cded6"]
            if sort_by == 'name':
                ax.set_title(self.region + ' amino acids', size=20)
                for i, amino_acid in enumerate(sorted(amino_acid_index.keys())):
                    c = colors[i]
                    ax.bar(position, aa[amino_acid_index[amino_acid], position], bottom=previous,
                           label=amino_acid, color=c)
                    previous += aa[amino_acid_index[amino_acid], position]
                lgd = ax.legend(sorted(amino_acid_index.keys()), loc='center left', bbox_to_anchor=(1, 0.5),
                                prop={"size": 16})

            elif sort_by == 'hydropathy':
                colors = ['b', 'r', 'k']
                ax.set_title(self.region + ' amino acid hydropathy', size=20)
                for i, prop_i in enumerate(['Hydrophilic', 'Moderate', 'Hydrophobic']):
                    c = colors[i]
                    ax.bar(position, hyd[i, position], bottom=previous, label=prop_i, color=c)
                    previous += hyd[i, position]
                lgd = ax.legend(['Hydrophilic', 'Moderate', 'Hydrophobic'], loc='center left', bbox_to_anchor=(1, 0.5),
                                prop={"size": 16})

            else:
                colors = ['b', 'r', 'k']
                ax.set_title(self.region + ' amino acid charge', size=20)
                for i, prop_i in enumerate(['Negative', 'Positive', 'Neutral']):
                    c = colors[i]
                    ax.bar(position, chg[i, position], bottom=previous, label=prop_i, color=c)
                    previous += chg[i, position]
                lgd = ax.legend(['Negative', 'Positive', 'Neutral'], loc='center left', bbox_to_anchor=(1, 0.5),
                                prop={"size": 16})

        if display_count:

            if self._aa_count.sum(0).max() < 99:

                # so that the text is always at the same distance of the bar
                if aa.max() > 1:
                    shift = aa.max() * 0.05
                else:
                    shift = 0.02

                for position in range(aa.shape[1]):
                    if self._aa_count[:, position].sum() > 9:
                        ax.text(x=position, y=aa[:, position].sum()+2*shift,
                                s=str(int(self._aa_count[:, position].sum())), rotation=45)
                    else:
                        ax.text(x=position, y=aa[:, position].sum()+shift,
                                s=str(int(self._aa_count[:, position].sum())), rotation=45)

            else:

                if aa.max() > 1:
                    shift = aa.max() * 0.075
                else:
                    shift = 0.02

                for position in range(aa.shape[1]):
                    if self._aa_count[:, position].sum() < 99:
                        ax.text(x=position-0.2, y=aa[:, position].sum() + shift, rotation=45,
                                s=str(int(self._aa_count[:, position].sum())))
                    elif 99 < self._aa_count[:, position].sum() < 999:
                        ax.text(x=position-0.2, y=aa[:, position].sum() + 2 * shift, rotation=45,
                                s=str(int(self._aa_count[:, position].sum())))
                    else:
                        ax.text(x=position-0.2, y=aa[:, position].sum() + 4 * shift, rotation=45,
                                s=str(int(self._aa_count[:, position].sum())))

        if normalize:
            ax.set_ylabel('Frequency', size=16)
        else:
            ax.set_ylabel('Count', size=16)

        ax.set_xticks(np.arange(len(self._numbering)) + 0.3)
        ax.set_xticklabels(self._numbering, rotation=60)
        ax.set_xlabel('Position', size=16)
        ax.set_ylim([0, aa.sum(0).max()*1.1])
        ax.margins(0.02)
        ax.grid(axis='x')

        ipython_config = PythonConfig()
        ipython_config.get_ipython_info()
        if ipython_config.backend == 'notebook' and notebook_plot:
            ax.plot()
        else:
            fig.savefig(path.join(plot_path, plot_name), bbox_extra_artists=(lgd,), bbox_inches='tight')
            plt.close(fig)

import re

import numpy as np

from src.utils.data_loader import DataLoader
from src.utils.downloads import Download


class Antibody:
    """
    TODO: write description
    """

    def __init__(self, sequence='', name='', numbering=[]):
        self.raw_sequence = sequence
        self.sequence = sequence.replace('-', '').upper()
        self.name = name
        self.numbering = numbering
        self.hydrophobicity_matrix = np.array([])
        self.chain = ''
        self.mw = 0

    def load(self):
        """
        Generates all the data available:
        - Antibody Numbering
        - Hydrophobicity matrix
        :return:

        """
        self.numbering, self.chain = self.ab_numbering()
        self.hydrophobicity_matrix = self.ab_hydrophobicity_matrix()
        self.mw = self.ab_molecular_weight()

    def ab_numbering(self, server='abysis', numbering_scheme='chothia'):
        # type: (str, str) -> object

        available_numbering_schemes = ['chothia', 'chothia_ext', 'kabath']
        available_servers = ['abysis']

        assert numbering_scheme.lower() in available_numbering_schemes, \
            "Unknown Numbering scheme: {}. \
            Numbering schemes available: {}".format(numbering_scheme,
                                                    ', '.join(available_numbering_schemes))

        assert server in available_servers, "Unknown server: {}. \
            Available servers: {}".format(server, ' ,'.join(available_servers))

        numbering = get_ab_numbering(self.sequence, server, numbering_scheme)
        chain = ''

        if numbering == ['']:
            print('Could not apply numbering scheme on provided sequence')
            return 'NA', 'NA'

        elif numbering[0][0] == 'H':
                chain = 'Heavy'
        elif numbering[0][0] == 'L':
                chain = 'Light'

        return numbering, chain

    def ab_hydrophobicity_matrix(self, hydrophobicity_scores='ew'):

        # check if all the required parameters are in order
        available_hydrophobicity_scores = ['kd', 'ww', 'hh', 'mf', 'ew']

        if isinstance(hydrophobicity_scores, str):
            assert hydrophobicity_scores in available_hydrophobicity_scores, \
                "Chosen hydrophobicity scores ({}) not available. \
                Available hydrophobicity scores: {}".format(
                    hydrophobicity_scores, ' ,'.join(available_hydrophobicity_scores)
                )
        if self.chain == 'Light':
            data_loader = DataLoader(numbering='LightChothiaWithCDR')
            whole_sequence = data_loader.get_data()

        elif self.chain == 'Heavy':
            data_loader = DataLoader(numbering='HeavyChothiaWithCDR')
            whole_sequence = data_loader.get_data()

        else:
            self.hydrophobicity_matrix = 'NA'
            print('Could not calculate the hydrophobicity matrix of the \
                  the provided sequence')
            return np.array([])

        # get the dictionary with the hydrophobicity scores
        data_loader = DataLoader(amino_acid_property=['hydrophobicity', hydrophobicity_scores + 'Hydrophobicity'])
        aa_hydrophobicity_scores = data_loader.get_data()

        return calculate_hydrophobicity_matrix(whole_sequence=whole_sequence, numbering=self.numbering,
                                               aa_hydrophobicity_scores=aa_hydrophobicity_scores,
                                               sequence=self.sequence)

    def ab_molecular_weight(self):

        data_loader = DataLoader(amino_acid_property=['MolecularWeight', 'average'])
        mw_dict = data_loader.get_data()

        return calculate_mw(self.sequence, mw_dict)


def get_ab_numbering(sequence, server, numbering_scheme):

    """

    :rtype: list
    """
    # check which server to use to get numbering
    if server.lower() == 'abysis':
        # find out which numbering scheme to use
        if numbering_scheme.lower() == 'chothia':
            scheme = '-c'
        elif numbering_scheme.lower() == 'chotia_ext':
            scheme = '-a'
        else:
            scheme = '-k'

        url = 'http://www.bioinf.org.uk/cgi-bin/abnum/abnum.pl?plain=1&aaseq={}&scheme={}'.format(sequence,
                                                                                                  scheme)
        numbering_table = Download(url, verbose=False)
        numbering_table.download()

        parsed_numbering_table = re.findall("[\S| ]+", numbering_table.html)

        numbering = [x[:-2] for x in parsed_numbering_table if x[-1] != '-']

        # TODO: add more server options
    else:
        numbering = ['']

    return numbering


def calculate_hydrophobicity_matrix(whole_sequence, numbering, aa_hydrophobicity_scores, sequence):

    # instantiate numpy array
    hydrophobicity_matrix = np.zeros(len(whole_sequence))

    for i, position in enumerate(whole_sequence):

        if position not in numbering:
            hydrophobicity_matrix[i] = 0

        else:
            position_in_data = numbering.index(position)
            hydrophobicity_matrix[i] = aa_hydrophobicity_scores[sequence[position_in_data]]

    return hydrophobicity_matrix


def calculate_mw(sequence, mw_dict):

    return sum(mw_dict[x] for x in sequence) - (len(sequence) - 1) * mw_dict['water']

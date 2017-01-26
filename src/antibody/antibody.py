import re

import numpy as np

from src.utils.data_loader import DataLoader
from src.utils.downloads import Download


class Antibody:
    """
    TODO: write description
    """

    def __init__(self, sequence='', name='', numbering=None):
        self.raw_sequence = sequence
        self.sequence = sequence.replace('-', '').upper()
        self.name = name
        self.numbering = numbering
        self.hydrophobicity_matrix = np.array([])
        self.chain = ''
        self.mw = 0
        self.pI = 0

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
        self.pI = self.ab_pi()

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

    def ab_molecular_weight(self, monoisotopic=False):

        if monoisotopic:
            data_loader = DataLoader(amino_acid_property=['MolecularWeight', 'average'])
        else:
            data_loader = DataLoader(amino_acid_property=['MolecularWeight', 'monoisotopic'])
        mw_dict = data_loader.get_data()

        return calculate_mw(self.sequence, mw_dict)

    def ab_pi(self, pi_database='Wikipedia'):

        available_pi_databases = ["EMBOSS", "DTASetect", "Solomon", "Sillero", "Rodwell", "Wikipedia", "Lehninger",
                                  "Grimsley"]
        assert pi_database in available_pi_databases, \
            "Selected pI database {} not available. Available databases: {}".format(pi_database,
                                                                                    ' ,'.join(available_pi_databases))

        data_loader = DataLoader(amino_acid_property=['pI', pi_database])
        pi_data = data_loader.get_data()

        return calculate_pi(sequence=self.sequence, pi_data=pi_data)


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


def calculate_pi(sequence, pi_data):

    # algorithm from http://isoelectric.ovh.org/files/practise-isoelectric-point.html

    # count number of D, E, C, Y, H, K, R
    d_count = sequence.count('D')
    e_count = sequence.count('E')
    c_count = sequence.count('C')
    y_count = sequence.count('Y')
    h_count = sequence.count('H')
    k_count = sequence.count('K')
    r_count = sequence.count('R')

    # initiate value of pH and nq (any number above 0)
    nq = 10
    pH = 0
    # define precision
    delta = 0.01

    while nq > 0:

        if pH >= 14:
            print('Could not calculate pI')
            break

        # qn1, qn2, qn3, qn4, qn5, qp1, qp2, qp3, qp4
        qn1 = -1 / (1 + 10 ** (pi_data['COOH'] - pH))  # C-terminus charge
        qn2 = - d_count / (1 + 10 ** (pi_data['D'] - pH))  # D charge
        qn3 = - e_count / (1 + 10 ** (pi_data['E'] - pH))  # E charge
        qn4 = - c_count / (1 + 10 ** (pi_data['C'] - pH))  # C charge
        qn5 = - y_count / (1 + 10 ** (pi_data['Y'] - pH))  # Y charge
        qp1 = h_count / (1 + 10 ** (pH - pi_data['H']))  # H charge
        qp2 = 1 / (1 + 10 ** (pH - pi_data['NH2']))  # N-terminus charge
        qp3 = k_count / (1 + 10 ** (pH - pi_data['K']))  # K charge
        qp4 = r_count / (1 + 10 ** (pH - pi_data['R']))  # R charge

        nq = qn1 + qn2 + qn3 + qn4 + qn5 + qp1 + qp2 + qp3 + qp4

        # update pH
        pH += delta

    return pH

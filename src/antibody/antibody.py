import re
import numpy as np
from ..utils import DataLoader, Download
import logging

# setting up debugging messages
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)


class Antibody:
    """
    TODO: write description
    """

    def __init__(self, sequence='', name='', numbering=None):
        self._raw_sequence = sequence.upper()
        self.sequence = self._raw_sequence.replace('-', '')
        self.name = name
        self.numbering = numbering
        self.hydrophobicity_matrix = np.array([])
        self.chain = ''
        self.mw = 0
        self.pI = 0
        self.cdr = [0, 0, 0]
        self._numbering_scheme = 'chothia'

    def load(self):
        """
        Generates all the data:
        - Antibody Numbering
        - Hydrophobicity matrix
        - Molecular weight
        - pI

        All the data is then stored in its respective attributes

        :return:

        """
        try:
            self.numbering, self.chain = self.ab_numbering()
            self.hydrophobicity_matrix = self.ab_hydrophobicity_matrix()
            self.mw = self.ab_molecular_weight()
            self.pI = self.ab_pi()
            self.cdr = self.ab_cdr()
        except ValueError:
            self.numbering = 'NA'
            self.chain = 'NA'
            self.hydrophobicity_matrix = 'NA'
            self.mw = 'NA'
            self.pI = 'NA'
            self.cdr = 'NA'

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

        # store the numbering scheme used for reference in other methods
        self._numbering_scheme = numbering_scheme

        numbering = get_ab_numbering(self.sequence, server, numbering_scheme)
        chain = ''

        if numbering == ['']:
            print('Could not apply numbering scheme on provided sequence')
            return 'NA', 'NA'

        elif numbering[0][0] == 'H':
                chain = 'heavy'
        elif numbering[0][0] == 'L':
                chain = 'light'

        return numbering, chain

    def ab_numbering_as_pandas(self, server='abysis', numbering_scheme='chothia'):
        pass

    def ab_hydrophobicity_matrix(self, hydrophobicity_scores='ew', include_cdr=True):

        # check if all the required parameters are in order
        available_hydrophobicity_scores = ['kd', 'ww', 'hh', 'mf', 'ew']

        if isinstance(hydrophobicity_scores, str):
            assert hydrophobicity_scores in available_hydrophobicity_scores, \
                "Chosen hydrophobicity scores ({}) not available. \
                Available hydrophobicity scores: {}".format(
                    hydrophobicity_scores, ' ,'.join(available_hydrophobicity_scores)
                )

        if self.chain == '':
            self.chain, self.numbering = self.ab_numbering()
        if self.chain == 'NA':
            raise ValueError("Could not determine chain type")

        data_loader = DataLoader(data_type='NumberingSchemes',
                                 data=[self._numbering_scheme, self.chain])
        whole_sequence_dict = data_loader.get_data()

        if include_cdr:
            whole_sequence = whole_sequence_dict['withCDR']
        else:
            whole_sequence = whole_sequence_dict['noCDR']

        # get the dictionary with the hydrophobicity scores
        data_loader = DataLoader(data_type='AminoAcidProperties',
                                 data=['hydrophobicity', hydrophobicity_scores + 'Hydrophobicity'])
        aa_hydrophobicity_scores = data_loader.get_data()

        return calculate_hydrophobicity_matrix(whole_sequence=whole_sequence, numbering=self.numbering,
                                               aa_hydrophobicity_scores=aa_hydrophobicity_scores,
                                               sequence=self.sequence)

    def ab_cdr(self):

        if self.numbering is None:
            self.numbering, self.chain = self.ab_numbering()

        if self.numbering == 'NA':
            raise ValueError("Cannot return CDR length without the antibody numbering information")

        data_loader = DataLoader(data_type='CDR_positions', data=[self._numbering_scheme, self.chain])
        cdr_positions = data_loader.get_data()

        return calculate_cdr(numbering=self.numbering, cdr_positions=cdr_positions)

    def ab_molecular_weight(self, monoisotopic=False):

        if monoisotopic:
            data_loader = DataLoader(data_type='AminoAcidProperties',
                                     data=['MolecularWeight', 'average'])
        else:
            data_loader = DataLoader(data_type='AminoAcidProperties',
                                     data=['MolecularWeight', 'monoisotopic'])
        mw_dict = data_loader.get_data()

        return calculate_mw(self.sequence, mw_dict)

    def ab_pi(self, pi_database='Wikipedia'):

        available_pi_databases = ["EMBOSS", "DTASetect", "Solomon", "Sillero", "Rodwell", "Wikipedia", "Lehninger",
                                  "Grimsley"]
        assert pi_database in available_pi_databases, \
            "Selected pI database {} not available. Available databases: {}".format(pi_database,
                                                                                    ' ,'.join(available_pi_databases))

        data_loader = DataLoader(data_type='AminoAcidProperties',
                                 data=['pI', pi_database])
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

        if numbering_table.html.replace("\n", '') == 'Warning: Unable to number sequence' or len(
                numbering_table.html.replace("\n", '')) == 0:
            raise ValueError("Unable to number sequence")

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

    # algorithm implemented from http://isoelectric.ovh.org/files/practise-isoelectric-point.html

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
    ph = 0
    # define precision
    delta = 0.01

    while nq > 0:

        if ph >= 14:
            raise Exception("Could not calculate pI (pH reached above 14)")

        # qn1, qn2, qn3, qn4, qn5, qp1, qp2, qp3, qp4
        qn1 = -1 / (1 + 10 ** (pi_data['COOH'] - ph))  # C-terminus charge
        qn2 = - d_count / (1 + 10 ** (pi_data['D'] - ph))  # D charge
        qn3 = - e_count / (1 + 10 ** (pi_data['E'] - ph))  # E charge
        qn4 = - c_count / (1 + 10 ** (pi_data['C'] - ph))  # C charge
        qn5 = - y_count / (1 + 10 ** (pi_data['Y'] - ph))  # Y charge
        qp1 = h_count / (1 + 10 ** (ph - pi_data['H']))  # H charge
        qp2 = 1 / (1 + 10 ** (ph - pi_data['NH2']))  # N-terminus charge
        qp3 = k_count / (1 + 10 ** (ph - pi_data['K']))  # K charge
        qp4 = r_count / (1 + 10 ** (ph - pi_data['R']))  # R charge

        nq = qn1 + qn2 + qn3 + qn4 + qn5 + qp1 + qp2 + qp3 + qp4

        # update pH
        ph += delta

    return ph


def calculate_cdr(numbering, cdr_positions):

    cdr_values = list()

    cdrs = ['CDR1', 'CDR2', 'CDR3']

    for i, cdr in enumerate(cdrs):

        cdr_positions_i = cdr_positions[cdr]
        cdr_value_i = 0

        for position in numbering:
            if position in cdr_positions_i:
                cdr_value_i += 1

        cdr_values.append(cdr_value_i)

    return cdr_values

import re
import numpy as np
from ..utils import DataLoader, Download
import logging
import pandas as pd
from .helper_functions import numbering_table_sequences, numbering_table_region, numbering_table_multiindex

# setting up debugging messages
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

available_numbering_schemes = ['chothia', 'chothia_ext', 'kabath']
available_servers = ['abysis']
available_hydrophobicity_scores = ['kd', 'ww', 'hh', 'mf', 'ew']
available_pi_databases = ["EMBOSS", "DTASetect", "Solomon", "Sillero", "Rodwell", "Wikipedia", "Lehninger",
                          "Grimsley"]


class Chain:
    # TODO: write description
    """
    """

    def __init__(self, sequence='', name='Chain1', numbering=None, numbering_scheme='chothia'):
        self._raw_sequence = sequence.upper()
        self._sequence = self._raw_sequence.replace('-', '')
        self._name = name
        self.numbering = numbering
        self.hydrophobicity_matrix = np.array([])
        self._chain = ''
        self.mw = 0
        self.pI = 0
        self.cdr = [0, 0, 0]
        self._numbering_scheme = numbering_scheme
        self._loading_status = 'Not Loaded'
        self.germline_identity = dict()
        self.germline = tuple()

    def load(self):
        """
        Generates all the data:
        - Chain Numbering
        - Hydrophobicity matrix
        - Molecular weight
        - pI

        All the data is then stored in its respective attributes

        :return:

        """
        if self._loading_status in ['Failed', 'Not Loaded']:
            try:
                self.numbering = self.ab_numbering()
                self._loading_status = 'Loaded'
                self.load()

            except ValueError:
                self._loading_status = 'Failed'

        elif self._loading_status == 'Loaded':
            self.hydrophobicity_matrix = self.ab_hydrophobicity_matrix()
            self.mw = self.ab_molecular_weight()
            self.pI = self.ab_pi()
            self.cdr = self.ab_regions()

        else:
            # this should never happen...
            raise ValueError("Unknown loading status")

    def ab_numbering(self, server='abysis', numbering_scheme='chothia'):

        assert numbering_scheme.lower() in available_numbering_schemes, \
            "Unknown Numbering scheme: {}. \
            Numbering schemes available: {}".format(numbering_scheme,
                                                    ', '.join(available_numbering_schemes))

        assert server in available_servers, "Unknown server: {}. \
            Available servers: {}".format(server, ' ,'.join(available_servers))

        # store the numbering scheme used for reference in other methods
        self._numbering_scheme = numbering_scheme

        # store the amino positions/numbering in a list -> len(numbering) == len(self._sequence)
        numbering = get_ab_numbering(self._sequence, server, numbering_scheme)

        if numbering == ['']:
            self._loading_status = 'Failed'
            return 'NA'

        elif numbering[0][0] == 'H':
            self._chain = 'heavy'
        elif numbering[0][0] == 'L':
            self._chain = 'light'

        return numbering

    def ab_numbering_table(self, as_array=False, replacement='-', region='all'):

        """

        :param region: 
        :param as_array: if True returns numpy.array object, if False returns a pandas.DataFrame
        :param replacement: value to replace empty positions
        :return:
        """

        region = numbering_table_region(region=region)

        # if the object has not been loaded successfully yet need to try and get the numbering scheme using
        # ab_numbering method
        if self._loading_status in ['Not Loaded', 'Failed']:
            self.numbering = self.ab_numbering()

        whole_sequence_dict, whole_sequence = numbering_table_sequences(region=region,
                                                                        numbering_scheme=self._numbering_scheme,
                                                                        chain=self._chain)

        # now that all the prep has been done we can extract the position from each amino acid
        # according to the numbering scheme
        data = np.empty((len(whole_sequence)), dtype=str)
        for i, position in enumerate(whole_sequence):
            if position in self.numbering:
                data[i] = self._sequence[self.numbering.index(position)]
            else:
                # if there is no amino acid in the sequence that corresponds to position i we just replace it by
                # the replacement value, which is by default '-'
                data[i] = replacement

        if as_array:
            # return the data as numpy.array -> much faster than creating a pandas.DataFrame
            return data
        else:

            # return the data as a pandas.DataFrame -> it's slower but looks nicer and makes it easier to get
            # the data of interest

            multi_index = numbering_table_multiindex(region=region,
                                                     whole_sequence_dict=whole_sequence_dict)

            # create the DataFrame and assign the columns and index names
            data = pd.DataFrame(data=data).T
            data.columns = multi_index
            data.index = [self._name]
            return data

    def ab_hydrophobicity_matrix(self, hydrophobicity_scores='ew'):

        # check if all the required parameters are in order
        if isinstance(hydrophobicity_scores, str):
            assert hydrophobicity_scores in available_hydrophobicity_scores, \
                "Chosen hydrophobicity scores ({}) not available. \
                Available hydrophobicity scores: {}".format(
                    hydrophobicity_scores, ' ,'.join(available_hydrophobicity_scores)
                )

        if self._chain == '':
            self._chain, self.numbering = self.ab_numbering()
        if self._chain == 'NA':
            raise ValueError("Could not determine chain type")

        data_loader = DataLoader(data_type='NumberingSchemes',
                                 data=[self._numbering_scheme, self._chain])
        whole_sequence_dict = data_loader.get_data()

        # whole_sequence is a list with all the amino acid positions in the selected numbering scheme
        whole_sequence = whole_sequence_dict

        # get the dictionary with the hydrophobicity scores
        data_loader = DataLoader(data_type='AminoAcidProperties',
                                 data=['hydrophobicity', hydrophobicity_scores + 'Hydrophobicity'])
        aa_hydrophobicity_scores = data_loader.get_data()

        return calculate_hydrophobicity_matrix(whole_sequence=whole_sequence, numbering=self.numbering,
                                               aa_hydrophobicity_scores=aa_hydrophobicity_scores,
                                               sequence=self._sequence)

    def ab_regions(self):

        """
        method to determine Chain regions (CDR and Framework) of each amino acid in sequence

        :return:
        """

        if self.numbering is None:
            self.numbering, self._chain = self.ab_numbering()

        if self.numbering == 'NA':
            raise ValueError("Cannot return CDR positions without the antibody numbering information")

        data_loader = DataLoader(data_type='CDR_positions', data=[self._numbering_scheme, self._chain])
        cdr_positions = data_loader.get_data()

        data_loader = DataLoader(data_type='Framework_positions', data=[self._numbering_scheme, self._chain])
        framework_position = data_loader.get_data()

        return calculate_cdr(numbering=self.numbering, cdr_positions=cdr_positions,
                             framework_positions=framework_position)

    def ab_molecular_weight(self, monoisotopic=False):

        if monoisotopic:
            data_loader = DataLoader(data_type='AminoAcidProperties',
                                     data=['MolecularWeight', 'average'])
        else:
            data_loader = DataLoader(data_type='AminoAcidProperties',
                                     data=['MolecularWeight', 'monoisotopic'])
        mw_dict = data_loader.get_data()

        return calculate_mw(self._sequence, mw_dict)

    def ab_pi(self, pi_database='Wikipedia'):

        assert pi_database in available_pi_databases, \
            "Selected pI database {} not available. Available databases: {}".format(pi_database,
                                                                                    ' ,'.join(available_pi_databases))

        data_loader = DataLoader(data_type='AminoAcidProperties',
                                 data=['pI', pi_database])
        pi_data = data_loader.get_data()

        return calculate_pi(sequence=self._sequence, pi_data=pi_data)

    def ab_ec(self, extinction_coefficient_database='Standard', reduced=False, normalise=False, **kwargs):

        if reduced:
            extinction_coefficient_database += '_reduced'

        data_loader = DataLoader(data_type='AminoAcidProperties', data=['ExtinctionCoefficient',
                                                                        extinction_coefficient_database])

        ec_data = data_loader.get_data()

        if normalise:
            return calculate_ec(sequence=self._sequence, ec_data=ec_data) / self.ab_molecular_weight(**kwargs)
        else:
            return calculate_ec(sequence=self._sequence, ec_data=ec_data)

    def ab_format(self):
        return {"name": self._name, "sequence": self._sequence, "numbering": self.numbering, "chain": self._chain,
                "MW": self.mw, "CDR": self.cdr, "numbering_scheme": self._numbering_scheme, "pI": self.pI}

    def ab_charge(self, align=True, ph=7.4, pka_database='Wikipedia'):

        """
        Method to calculate the charges for each amino acid of antibody
        :param pka_database: 
        :param ph: 
        :param align: if set to True an alignment will be performed, 
                      if it hasn't been done already using the ab_numbering method
                        
        :return: array with amino acid charges
        """

        assert pka_database in available_pi_databases, \
            "Selected pI database {} not available. Available databases: {}".format(pka_database,
                                                                                    ' ,'.join(available_pi_databases))

        data_loader = DataLoader(data_type='AminoAcidProperties',
                                 data=['pI', pka_database])
        pka_data = data_loader.get_data()

        if align:
            sequence = self.ab_numbering_table(as_array=True)
        else:
            sequence = list(self.sequence)

        return np.array([amino_acid_charge(x, ph, pka_data) for x in sequence])

    def ab_total_charge(self, ph=7.4, pka_database='Wikipedia'):

        assert pka_database in available_pi_databases, \
            "Selected pI database {} not available. Available databases: {}".format(pka_database,
                                                                                    ' ,'.join(available_pi_databases))

        data_loader = DataLoader(data_type='AminoAcidProperties',
                                 data=['pI', pka_database])
        pka_data = data_loader.get_data()

        return calculate_charge(sequence=self._sequence, ph=ph, pka_values=pka_data)

    @property
    def chain(self):
        return self._chain

    @property
    def name(self):
        return self._name

    @property
    def sequence(self):
        return self._sequence

    @property
    def status(self):
        return self._loading_status

    @property
    def numbering_scheme(self):
        return self._numbering_scheme

    @property
    def insertions(self):
        return [(number, seq) for seq, number in zip(self.sequence, self.numbering) if number[-1].isdigit() is False]

    def _string_summary_basic(self):
        return "abpytools.Chain Name: {}, Chain type: {}, Sequence length: {}, Status: {}".format(self._name,
                                                                                                  self._chain,
                                                                                                  len(self._sequence),
                                                                                                  self._loading_status)

    def __repr__(self):
        return "<%s at 0x%02x>" % (self._string_summary_basic(), id(self))

    def __len__(self):
        return len(self.sequence)


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

        # prepare the url string to query server
        url = 'http://www.bioinf.org.uk/cgi-bin/abnum/abnum.pl?plain=1&aaseq={}&scheme={}'.format(sequence,
                                                                                                  scheme)
        # use the Download class from utils to get output
        numbering_table = Download(url, verbose=False)
        try:
            numbering_table.download()
        except ValueError:
            raise ValueError("Check the internet connection.")

        # check whether the server returned an error
        if numbering_table.html.replace("\n", '') == 'Warning: Unable to number sequence' or len(
                numbering_table.html.replace("\n", '')) == 0:
            raise ValueError("Unable to number sequence")

        # parse the results
        parsed_numbering_table = re.findall("[\S| ]+", numbering_table.html)

        # get the numbering from the parsed table
        numbering = [x[:-2] for x in parsed_numbering_table if x[-1] != '-']

        # TODO: add more server options
    else:
        numbering = ['']

    return numbering


def calculate_hydrophobicity_matrix(whole_sequence, numbering, aa_hydrophobicity_scores, sequence):
    # instantiate numpy array (whole sequence includes all the amino acid positions of the VH/VL, even the ones
    # that aren't occupied -> these will be filled with zeros
    # hydrophobicity_matrix = np.zeros(len(whole_sequence))
    #
    # # iterate through each position
    # for i, position in enumerate(whole_sequence):
    #
    #     if position in numbering:
    #         position_in_data = numbering.index(position)
    #         hydrophobicity_matrix[i] = aa_hydrophobicity_scores[sequence[position_in_data]]
    #
    # return hydrophobicity_matrix
    # same thing as above but in a list comprehension
    return np.array([aa_hydrophobicity_scores[sequence[numbering.index(x)]] if x in numbering
                     else 0 for x in whole_sequence])


def calculate_mw(sequence, mw_data):
    return sum(mw_data[x] for x in sequence) - (len(sequence) - 1) * mw_data['water']


def calculate_ec(sequence, ec_data):
    # ϵ280 = nW x 5,500 + nY x 1,490 + nC x 125
    n_W = sequence.count('W')
    n_Y = sequence.count('Y')
    n_C = sequence.count('C')
    return n_W * ec_data['W'] + n_Y * ec_data['Y'] + n_C * ec_data['C']


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


def calculate_cdr(numbering, cdr_positions, framework_positions):
    """

    :param numbering:
    :param cdr_positions:
    :param framework_positions:
    :return:
    """

    cdrs = {'CDR1': list(),
            'CDR2': list(),
            'CDR3': list()}

    frameworks = {'FR1': list(),
                  'FR2': list(),
                  'FR3': list(),
                  'FR4': list()}

    for cdr in cdrs.keys():

        cdr_positions_i = cdr_positions[cdr]

        for i, position in enumerate(numbering):
            if position in cdr_positions_i:
                cdrs[cdr].append(i)

    for framework in frameworks.keys():

        framework_position_i = framework_positions[framework]

        for i, position in enumerate(numbering):
            if position in framework_position_i:
                frameworks[framework].append(i)

    return cdrs, frameworks


def amino_acid_charge(amino_acid, ph, pka_values):
    if amino_acid in ['D', 'E', 'C', 'Y']:
        return -1 / (1 + 10 ** (pka_values[amino_acid] - ph))
    elif amino_acid in ['K', 'R', 'H']:
        return 1 / (1 + 10 ** (ph - pka_values[amino_acid]))
    else:
        return 0


def calculate_charge(sequence, ph, pka_values):
    # This calculation would make more sense but is slower (~1.5-2x)
    # cooh = -1 / (1 + 10 ** (pka_values['COOH'] - ph))
    # nh2 = 1 / (1 + 10 ** (ph - pka_values['NH2']))
    #
    # return sum([amino_acid_charge(x, ph, pka_values) for x in list(sequence)]) + cooh + nh2

    # Faster implementation
    # count number of D, E, C, Y, H, K, R
    d_count = sequence.count('D')
    e_count = sequence.count('E')
    c_count = sequence.count('C')
    y_count = sequence.count('Y')
    h_count = sequence.count('H')
    k_count = sequence.count('K')
    r_count = sequence.count('R')

    # qn1, qn2, qn3, qn4, qn5, qp1, qp2, qp3, qp4
    qn1 = -1 / (1 + 10 ** (pka_values['COOH'] - ph))  # C-terminus charge
    qn2 = - d_count / (1 + 10 ** (pka_values['D'] - ph))  # D charge
    qn3 = - e_count / (1 + 10 ** (pka_values['E'] - ph))  # E charge
    qn4 = - c_count / (1 + 10 ** (pka_values['C'] - ph))  # C charge
    qn5 = - y_count / (1 + 10 ** (pka_values['Y'] - ph))  # Y charge
    qp1 = h_count / (1 + 10 ** (ph - pka_values['H']))  # H charge
    qp2 = 1 / (1 + 10 ** (ph - pka_values['NH2']))  # N-terminus charge
    qp3 = k_count / (1 + 10 ** (ph - pka_values['K']))  # K charge
    qp4 = r_count / (1 + 10 ** (ph - pka_values['R']))  # R charge

    nq = qn1 + qn2 + qn3 + qn4 + qn5 + qp1 + qp2 + qp3 + qp4

    return nq

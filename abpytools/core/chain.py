import re
import numpy as np
from ..utils import DataLoader, Download, NumberingException
import pandas as pd
from .helper_functions import numbering_table_sequences, numbering_table_region, numbering_table_multiindex
from . import Cache
from .flags import *


class Chain:
    """The Chain object represent a single chain variable fragment (scFv) antibody.

    A scFv can be part of either the heavy or light chain of an antibody.
    The nature of the chain is determined by querying the sequence
    to the Abnum server, and is implemented with the Chain.ab_numbering() method.

    Attributes:
        numbering (list): the name of each position occupied by amino acids in sequence
        mw (float): the cached molecular weight
        pI (float): the cached isoelectric point of the sequence
        cdr (tuple): tuple with two dictionaries for CDR and FR with the index of the amino acids
        in each region
        germline_identity (dict):
    """

    def __init__(self, sequence, name='Chain1', numbering_scheme=NUMBERING_FLAGS.CHOTHIA):
        """
        The default Chain object constructor, which required a string
        representing a scFv sequence.

        Args:
            sequence (str): Amino acid sequence
            name (name): Name of sequence
            numbering_scheme (str): numbering scheme name to perform alignment

        Examples:
            Instantiate a Chain object with the default constructor
            >>> from abpytools.core import Chain
            >>> from abpytools.core.flags import *
            >>> chain = Chain(sequence='MYSEQUENCE', name='my_seq', numbering_scheme=NUMBERING_FLAGS.CHOTHIA)
        """
        self._raw_sequence = sequence.upper()
        self._sequence = self._raw_sequence.replace('-', '')
        self._aligned_sequence = None
        self._name = name
        self._chain = None
        self.numbering = None
        self.hydrophobicity_matrix = None
        self.mw = None
        self.pI = None
        self.cdr = None
        self._numbering_scheme = numbering_scheme
        self._loading_status = 'Not Loaded'
        self.germline_identity = dict()
        self.germline = tuple()
        self._cache = Cache(max_cache_size=10)

    @classmethod
    def load_from_string(cls, sequence, name='Chain1', numbering_scheme=NUMBERING_FLAGS.CHOTHIA):
        """
        Returns an instantiated Chain object from a sequence
        Args:
            sequence:
            name:
            numbering_scheme:

        Returns:

        """
        new_chain = cls(sequence=sequence, name=name, numbering_scheme=numbering_scheme)
        new_chain.load()
        return new_chain

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

        if self._loading_status in [NUMBERING_FLAGS.FAILED, NUMBERING_FLAGS.NOT_LOADED]:
            try:
                self.numbering = self.ab_numbering()
                self._loading_status = NUMBERING_FLAGS.LOADED
                self.load()

            except ValueError:
                self._loading_status = NUMBERING_FLAGS.FAILED

            except NumberingException:
                self._loading_status = NUMBERING_FLAGS.UNNUMBERED

        elif self._loading_status == NUMBERING_FLAGS.LOADED:
            self.hydrophobicity_matrix = self.ab_hydrophobicity_matrix()
            self.mw = self.ab_molecular_weight()
            self.pI = self.ab_pi()
            self.cdr = self.ab_regions()

        else:
            # this should never happen...
            raise ValueError("Unknown loading status")  # pragma: no cover

    @staticmethod
    def determine_chain_type(numbering):
        if numbering[0][0] == 'H':
            chain = CHAIN_FLAGS.HEAVY_CHAIN
        elif numbering[0][0] == 'L':
            chain = CHAIN_FLAGS.LIGHT_CHAIN
        else:
            # couldn't determine chain type
            chain = CHAIN_FLAGS.UNKNOWN_CHAIN  # pragma: no cover

        return chain

    def ab_numbering(self, server=OPTION_FLAGS.ABYSIS, **kwargs):
        """
        Return list

        Returns:
            list:
        """
        # store the amino positions/numbering in a list -> len(numbering) == len(self._sequence)
        numbering = get_ab_numbering(self._sequence, server, self._numbering_scheme, **kwargs)
        self._chain = self.determine_chain_type(numbering)

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
        if self._loading_status in [NUMBERING_FLAGS.NOT_LOADED,
                                    NUMBERING_FLAGS.FAILED]:
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

        self._aligned_sequence = data

        if as_array:
            # return the data as numpy.array -> much faster than creating a pandas.DataFrame
            return data.reshape(1, -1)
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

    def ab_hydrophobicity_matrix(self, hydrophobicity_scores=HYDROPHOBICITY_FLAGS.EW):

        # check if all the required parameters are in order
        if isinstance(hydrophobicity_scores, str):
            if hydrophobicity_scores not in OPTION_FLAGS.AVAILABLE_HYDROPHOBITY_SCORES:
                raise ValueError("Chosen hydrophobicity scores ({}) not available. \
                Available hydrophobicity scores: {}".format(
                    hydrophobicity_scores, ' ,'.join(OPTION_FLAGS.AVAILABLE_HYDROPHOBITY_SCORES)
                ))

        if self._loading_status == NUMBERING_FLAGS.NOT_LOADED:
            self.numbering = self.ab_numbering()
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

        return calculate_hydrophobicity_matrix(whole_sequence=whole_sequence,
                                               numbering=self.numbering,
                                               aa_hydrophobicity_scores=aa_hydrophobicity_scores,
                                               sequence=self._sequence)

    def ab_regions(self):

        """
        method to determine Chain regions (CDR and Framework) of each amino acid in sequence

        :return:
        """

        if 'cdrs' not in self._cache:

            if self._loading_status == NUMBERING_FLAGS.NOT_LOADED:
                self.numbering, self._chain = self.ab_numbering()

            if self.numbering == 'NA':
                raise ValueError("Cannot return CDR positions without the antibody numbering information")

            data_loader = DataLoader(data_type='CDR_positions', data=[self._numbering_scheme, self._chain])
            cdr_positions = data_loader.get_data()

            data_loader = DataLoader(data_type='Framework_positions', data=[self._numbering_scheme, self._chain])
            framework_position = data_loader.get_data()

            cdrs = calculate_cdr(numbering=self.numbering, cdr_positions=cdr_positions,
                                 framework_positions=framework_position)

            self._cache.update('cdrs', cdrs)

        return self._cache['cdrs']

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

        assert pi_database in OPTION_FLAGS.AVAILABLE_PI_VALUES, \
            "Selected pI database {} not available. " \
            "Available databases: {}".format(pi_database, ' ,'.join(
                OPTION_FLAGS.AVAILABLE_PI_VALUES))

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

        assert pka_database in OPTION_FLAGS.AVAILABLE_PI_VALUES, \
            "Selected pI database {} not available. " \
            "Available databases: {}".format(pka_database,
                                             ' ,'.join(OPTION_FLAGS.AVAILABLE_PI_VALUES))

        data_loader = DataLoader(data_type='AminoAcidProperties',
                                 data=['pI', pka_database])
        pka_data = data_loader.get_data()

        if align:
            # get the first (and only) row
            sequence = self.ab_numbering_table(as_array=True)[0]
        else:
            sequence = list(self.sequence)

        return np.array([amino_acid_charge(x, ph, pka_data) for x in sequence])

    def ab_total_charge(self, ph=7.4, pka_database=PI_FLAGS.WIKIPEDIA):

        assert pka_database in OPTION_FLAGS.AVAILABLE_PI_VALUES, \
            "Selected pI database {} not available. " \
            "Available databases: {}".format(pka_database,
                                             ' ,'.join(OPTION_FLAGS.AVAILABLE_PI_VALUES))

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

    def set_name(self, name):
        self._name = name

    @property
    def sequence(self):
        return self._sequence

    @property
    def aligned_sequence(self):
        if self._aligned_sequence is None:
            _ = self.ab_numbering_table(as_array=True, replacement='-')
        return self._aligned_sequence.tolist()

    @property
    def status(self):
        return self._loading_status

    @property
    def numbering_scheme(self):
        return self._numbering_scheme

    def _string_summary_basic(self):
        return "abpytools.Chain Name: {}, Chain type: {}, Sequence length: {}, Status: {}".format(
            self._name, self._chain, len(self._sequence),
            self._loading_status)  # pragma: no cover

    def __repr__(self):
        return "<%s at 0x%02x>" % (self._string_summary_basic(), id(self))  # pragma: no cover

    def __len__(self):
        return len(self.sequence)


def get_ab_numbering(sequence, server, numbering_scheme, timeout=30):
    """

    :rtype: list
    """
    # check which server to use to get numbering
    if server == OPTION_FLAGS.ABYSIS:
        # find out which numbering scheme to use
        if numbering_scheme == NUMBERING_FLAGS.CHOTHIA:
            scheme = '-c'
        elif numbering_scheme in (NUMBERING_FLAGS.CHOTHIA_EXT, NUMBERING_FLAGS.MARTIN):
            scheme = '-a'
        elif numbering_scheme == NUMBERING_FLAGS.KABAT:
            scheme = '-k'
        else:
            raise ValueError("{} numbering scheme is unknown.".format(numbering_scheme.capitalize()))

        # prepare the url string to query server
        url = f"http://www.bioinf.org.uk/abs/abnum/abnum.cgi?plain=1&aaseq={sequence}&scheme={scheme}"
        # use the Download class from utils to get output
        numbering_table = Download(url, verbose=False, timeout=timeout)
        try:
            numbering_table.download()
        except ValueError:
            raise ValueError("Check the internet connection.")

        # check whether the server returned an error
        if numbering_table.html.replace("\n", '') == 'Warning: Unable to number sequence' or len(
                numbering_table.html.replace("\n", '')) == 0:
            raise NumberingException("Unable to number sequence")

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
    # same thing as above but in a comprehension list
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

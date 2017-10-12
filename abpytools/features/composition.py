from collections import Counter, defaultdict
from itertools import  product
import re


aa_order = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

# this classification is based on the Shen paper and takes into account the aa side chain dipole and volume
# for more information check http://www.pnas.org/content/104/11/4337/suppl/DC1
aa_group = {'A': '0', 'G': '0', 'V': '0', 'I': '1', 'L': '1', 'F': '1', 'P': '1', 'Y': '2', 'M': '2', 'T': '2',
            'S': '2', 'H': '3', 'N': '3', 'Q': '3', 'W': '3', 'R': '4', 'K': '4', 'D': '5', 'E': '5', 'C': '6'}


def chou_pseudo_aa_composition(sequences):
    # M.K. Gupta , R. Niyogi & M. Misra (2013) An alignment-free method to find
    # similarity among protein sequences via the general form of Chou’s pseudo amino acid composition,
    # SAR and QSAR in Environmental Research, 24:7, 597-609,
    # DOI: 10.1080/1062936X.2013.773378

    # first the aa count
    aa_count_dict = [aa_composition(seq) for seq in sequences]

    # distance to first
    aa_distance_to_first_dict = [distance_to_first(x) for x in sequences]

    aa_distribution_dict = [aa_distribution(seq, aa_c, aa_dist) for seq, aa_c, aa_dist in
                            zip(sequences, aa_count_dict, aa_distance_to_first_dict)]

    # create lists with amino acids in the right order
    aa_count = [order_seq(aa_count_dict_i) for aa_count_dict_i in aa_count_dict]

    aa_distance_to_first = [order_seq(aa_distance_to_first_i) for aa_distance_to_first_i in aa_distance_to_first_dict]

    aa_dist = [order_seq(aa_distribution_dict_i) for aa_distribution_dict_i in aa_distribution_dict]

    return [x + y + z for x, y, z in zip(aa_count, aa_distance_to_first, aa_dist)]


def aa_composition(seq):
    return Counter(seq)


def aa_frequency(seq):
    aa_count = aa_composition(seq)
    total = sum(aa_count.values())
    return {key: value/total for key, value in aa_count.items()}


def distance_to_first(seq):
    return {x: sum([m.start() for m in re.finditer(x, seq)]) for x in aa_order}


def aa_distribution(seq, aa_count, aa_distance_to_first):
    aa_dist_dict = defaultdict(int)
    for i, aa in enumerate(seq):
        aa_dist_dict[aa] += (i - (aa_distance_to_first[aa] / aa_count[aa])) ** 2 / aa_count[aa]
    return aa_dist_dict


def order_seq(seq_dict):
    return [seq_dict[aa] if aa in seq_dict else 0 for aa in aa_order]


def triad_method(sequences):

    # as described in
    # Shen J. et al. (2006). Predicting protein–protein interactions based only on sequences information. PNAS,
    # 104(11), pp: 4337-4341.

    d_matrix = []

    for sequence in sequences:

        # start dictionary with all 343 triads/keys (7 classes ** 3)
        f_keys = [''.join(x) for x in product(['0', '1', '2', '3', '4', '5', '6'],
                                              ['0', '1', '2', '3', '4', '5', '6'],
                                              ['0', '1', '2', '3', '4', '5', '6'])]

        f_results = {x: 0 for x in f_keys}

        # classify aa in sequence
        v = [aa_group[aa] for aa in sequence]

        # get triads
        triads = [''.join(v[x:x + 3]) for x in range(len(v) - 2)]

        for triad in triads:
            f_results[triad] += 1

        # normalise values
        f_max = max(f_results.values())
        f_min = min(f_results.values())

        d = [(f_results[key] - f_min) / f_max for key in f_keys]

        # append 343 dimensional vector to d_matrix
        d_matrix.append(d)

    # d_matrix has shape (len(sequences), 343)
    return d_matrix


def side_chain_volume(sequences):
    pass


def auto_covariance(sequences):
    pass
from ..utils.math_utils import dot_product, magnitude
from math import acos
from ..core.cache import Cache
from .analysis_helper_functions import init_score_matrix


def cosine_distance(u, v):
    """
    returns the cosine distance between vectors u and v
    :param u:
    :param v:
    :return:
    """
    if u == v:
        return 0
    else:
        # this clips the values to be between 1 and -1, as there can be rounding errors
        return acos(max(min(dot_product(u, v) / (magnitude(u) * magnitude(v)), 1), -1))


def cosine_similarity(u, v):
    """
    returns the cosine similarity between vectors u and v
    :param u:
    :param v:
    :return:
    """
    return 1 - cosine_distance(u, v)


def hamming_distance(seq1, seq2):
    """
    returns the hamming distance between two sequences
    :param seq1:
    :param seq2:
    :return:
    """
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be equal length, instead got {} and {}".format(len(seq1), len(seq2)))
    return sum(aa1 != aa2 for aa1, aa2 in zip(seq1, seq2))


def levenshtein_distance(seq1, seq2):
    """
    levenshtein distance using dynamic programming
    """

    dist = init_score_matrix(seq_1=seq1, seq_2=seq2, indel=1)

    for col in range(1, len(dist[0])):
        for row in range(1, len(dist)):
            if seq2[row - 1] == seq1[col - 1]:
                cost = 0
            else:
                cost = 1
            dist[row][col] = min(dist[row - 1][col] + 1,  # deletion
                                 dist[row][col - 1] + 1,  # insertion
                                 dist[row - 1][col - 1] + cost)  # substitution

    return dist[-1][-1]

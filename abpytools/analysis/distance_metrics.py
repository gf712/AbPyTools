from .distance_metrics_ import cosine_distance_, hamming_distance_, levenshtein_distance_
from .analysis_helper_functions import init_score_matrix
# from math import acos
# from ..utils.math_utils import dot_product, magnitude


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
        # pure python: acos(max(min(dot_product(u, v) / (magnitude(u) * magnitude(v)), 1), -1))
        return cosine_distance_(u, v)


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
    # pure python: sum(aa1 != aa2 for aa1, aa2 in zip(seq1, seq2))
    return hamming_distance_(seq1, seq2)


def levenshtein_distance(seq1, seq2):
    """

    :param seq1:
    :param seq2:
    :return:
    """

    # pure python:
    # dist = init_score_matrix(seq_1=seq1, seq_2=seq2, indel=1)
    #
    # cols = len(dist[0])
    # rows = len(dist)
    #
    # for col in range(1, cols):
    #     for row in range(1, rows):
    #         if seq2[row - 1] == seq1[col - 1]:
    #             cost = 0
    #         else:
    #             cost = 1
    #         dist[row][col] = min(dist[row - 1][col] + 1,  # deletion
    #                              dist[row][col - 1] + 1,  # insertion
    #                              dist[row - 1][col - 1] + cost)  # substitution
    #
    # return dist[rows-1][cols-1]
    return levenshtein_distance_(seq1, seq2)

def euclidean_distance(u, v):
    """
    returns the euclidean distance
    :param u:
    :param v:
    :return:
    """
    return norm(u, v, degree=2)


def manhattan_distance(u, v):
    """
    returns the Manhattan distance
    :param u:
    :param v:
    :return:
    """
    return norm(u, v, degree=1)


def norm(u, v, degree=2):
    return sum([(abs(u_i - v_i))**degree for u_i, v_i in zip(u, v)]) ** (1/degree)

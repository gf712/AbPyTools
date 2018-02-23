import warnings
import _pickle as cPickle
from abpytools.home import Home
from ..utils.python_config import PythonConfig
import matplotlib.pyplot as plt

SUPPORTED_SUBSITUTION_MATRICES = ['BLOSUM45', 'BLOSUM62', 'BLOSUM80']


def load_alignment_algorithm(algorithm):
    if algorithm.lower() == 'needleman_wunsch':
        return needleman_wunsch

    else:
        raise ValueError("Unknown algorithm")


def load_substitution_matrix(substitution_matrix):

    abpytools_directory = Home().homedir

    if substitution_matrix in SUPPORTED_SUBSITUTION_MATRICES:
        with open('{}/data/{}.txt'.format(abpytools_directory, substitution_matrix), 'rb') as f:
            matrix = cPickle.load(f)
    else:
        raise ValueError("Unknown substitution matrix")

    return matrix


def needleman_wunsch(seq_1, seq_2, substitution_matrix, indel=-1):

    if indel >= 0:
        f = "Indel must be negative, setting indel to {}.".format(-indel)
        warnings.warn(f)
        indel = -indel

    # initialise matrix
    init_matrix = init_score_matrix(seq_1=seq_1, seq_2=seq_2,
                                    indel=indel)

    scores, traceback_matrix = calculate_scores(matrix=init_matrix,
                                                seq_1=seq_1,
                                                seq_2=seq_2,
                                                substitution_matrix=substitution_matrix,
                                                gap_penalty=indel)

    seq_2_aligned = traceback(traceback_matrix=traceback_matrix, seq_1=seq_1, seq_2=seq_2)

    return seq_2_aligned, scores[-1][-1]


def init_score_matrix(seq_1, seq_2, indel):
    """
    - score matrix initialisation with two sequences
    - pure python, i.e. no numpy
    Example init_score_matrix('SEND', 'AND', -1):
        [[0, -1, -2],
         [-1, 0, 0],
         [-2, 0, 0],
         [-3, 0, 0]]
    """

    init_matrix = [[x * indel] + [0] * len(seq_1) if x > 0 else list(range(0, (len(seq_1) + 1) * indel, indel)) for x in
                   range(len(seq_2) + 1)]

    return init_matrix


def calculate_scores(matrix, seq_1, seq_2, substitution_matrix, gap_penalty):
    traceback_matrix = [['up'] + [''] * len(seq_1) if x > 0 else ['left'] * (len(seq_1) + 1) for x in
                        range(len(seq_2) + 1)]
    traceback_matrix[0][0] = 'done'

    for i in range(1, len(matrix)):
        for j in range(1, len(matrix[i])):
            q_diag = matrix[i - 1][j - 1] + int(substitution_matrix[(seq_2[i - 1], seq_1[j - 1])])
            q_up = matrix[i - 1][j] + gap_penalty
            q_left = matrix[i][j - 1] + gap_penalty
            results = [q_diag, q_up, q_left]
            matrix[i][j] = max(results)
            max_index = results.index(matrix[i][j])
            if max_index == 0:
                traceback_matrix[i][j] = 'diag'
            elif max_index == 1:
                traceback_matrix[i][j] = 'up'
            else:
                traceback_matrix[i][j] = 'left'

    return matrix, traceback_matrix


def traceback(traceback_matrix, seq_1, seq_2):
    row = -1
    column = -1
    current = traceback_matrix[row][column]

    # iter_seq_1 = iter(seq_1[::-1])
    # iter_seq_2 = iter(seq_2[::-1])

    # seq_1 = seq_1[::-1]
    # seq_2 = seq_2[::-1]

    # aligned_seq_1 = list()
    aligned_seq_2 = list()

    while current != 'done':
        # is aligned
        if current == 'diag':
            # aligned_seq_1.append(next(iter_seq_1))
            # aligned_seq_2.append(next(iter_seq_2))
            # aligned_seq_1.append(seq_1[column])
            aligned_seq_2.append(seq_2[row])
            row -= 1
            column -= 1

        # leave a gap
        elif current == 'left':
            # aligned_seq_1.append(next(iter_seq_1))
            # aligned_seq_1.append(seq_1[column])
            aligned_seq_2.append('-')
            column -= 1

        # leave a gap
        else:
            # aligned_seq_1.append(next(iter_seq_1))
            # aligned_seq_1.append(seq_1[column])
            aligned_seq_2.append('-')
            row -= 1

        current = traceback_matrix[row][column]

    return ''.join(aligned_seq_2)[::-1]


def switch_interactive_mode(save=False):
    ipython_config = PythonConfig()
    if ipython_config.ipython_info == 'notebook' and save is False:
        if ipython_config.matplotlib_interactive is False:
            # turns on interactive mode
            plt.ion()

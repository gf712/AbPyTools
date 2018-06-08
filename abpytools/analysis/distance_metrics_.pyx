from abpytools.cython_extensions.convert_py_2_C cimport get_C_double_array_pointers, release_C_pointer, get_C_char_array_pointers
from abpytools.utils.math_utils_ cimport internal_vector_dot_product_, internal_vector_norm_
from libc.math cimport acos as acos_C
from libc.math cimport fmin as min_C
from libc.float cimport DBL_EPSILON
# from abpytools.utils.math_utils_ import Matrix2D
cimport numpy as np
import numpy as np


cpdef np.ndarray init_score_matrix_(str seq_1, str seq_2, int indel):
    """
    - score matrix initialisation with two sequences
    Example init_score_matrix('SEND', 'AND', -1):
        [[0, -1, -2],
         [-1, 0, 0],
         [-2, 0, 0],
         [-3, 0, 0]]
         
    Args:
        seq_1: 
        seq_2: 
        indel: 

    Returns:

    """

    init_matrix = np.array([[x * indel] + [0] * len(seq_1) if x > 0 else list(range(0, (len(seq_1) + 1) * indel, indel)) for x in
                            range(len(seq_2) + 1)])

    return init_matrix


cpdef double cosine_distance_(list u, list v):
    """
    
    Args:
        u: 
        v: 

    Returns:

    """
    cdef int size
    if len(u) == len(v):
        size = len(u)
    else:
        raise ValueError("Vector size mismatch")

    cdef double *u_ = get_C_double_array_pointers(u, size)
    cdef double *v_ = get_C_double_array_pointers(v, size)
    cdef double result
    cdef double lower
    cdef int norm = 2

    cdef double upper = internal_vector_dot_product_(u_, v_, size)

    if upper < DBL_EPSILON:
        # if the result is lower than the double rounding error just assume that it should be in fact 0
        result = 0

    else:
        lower = internal_vector_norm_(u_, norm, size) * internal_vector_norm_(v_, norm, size)
        result = acos_C(upper / lower)

    release_C_pointer(u_)
    release_C_pointer(v_)

    return result

cpdef int hamming_distance_(str seq1, str seq2):
    """
    
    Args:
        seq1: 
        seq2: 

    Returns:

    """

    cdef int size
    if len(seq1) == len(seq2):
        size = len(seq1)
    else:
        raise ValueError("Sequence size mismatch")

    cdef int i
    cdef int result = 0
    cdef char *seq1_ = get_C_char_array_pointers(seq1, size)
    cdef char *seq2_ = get_C_char_array_pointers(seq2, size)

    for i in range(size):
        if seq1_[i] != seq2_[i]:
            result+=1

    release_C_pointer(seq1_)
    release_C_pointer(seq2_)

    return result


cpdef levenshtein_distance(seq1, seq2):
    """
    
    Args:
        seq1: 
        seq2: 

    Returns:

    """
    cdef np.ndarray dist = init_score_matrix_(seq_1=seq1, seq_2=seq2, indel=1)
    cdef int cols = dist.shape[1]
    cdef int rows = dist.shape[0]
    cdef int col
    cdef int row
    cdef int cost

    for col in range(1, cols):
        for row in range(1, rows):
            if seq2[row - 1] == seq1[col - 1]:
                cost = 0
            else:
                cost = 1
            dist[row, col] = min_C(min_C(dist[row - 1, col] + 1,  # deletion
                                   dist[row, col - 1] + 1),  # insertion
                                   dist[row - 1, col - 1] + cost)  # substitution

    return dist[-1][-1]
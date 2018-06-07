from libc.math cimport pow as pow_C
from abpytools.cython_extensions.convert_py_2_C cimport get_C_double_array_pointers, release_C_pointer
import itertools


cdef class Matrix2D:

    def __init__(self, *rows):
        # check all rows have the same size
        self.rows = rows
        self.n_rows = len(rows)
        self.n_cols = len(rows[0])

        cdef int i

        for i in range(1, self.n_rows):
           if len(self.rows[i]) != self.n_cols:
               raise ValueError("All rows must have the same number of elements")

        # data represents elements in a contiguous array
        self.data = list(itertools.chain(*self.rows))
        # copy to a C array
        self.data_ = get_C_double_array_pointers(self.data, self.n_rows * self.n_cols)

    def __dealloc__(self):
        # releases C memory
        release_C_pointer(self.data_)

    def __getitem__(self, int index):

        return Matrix2D(self.data_[index])


# cdef Vector:



cdef double internal_vector_dot_product_(double *u, double *v, int size):
    """
    "Pure" C definition of dot product in Cython compiler to be only used in the backend
    
    Args:
        u: 
        v: 
        size: 

    Returns:

    """

    cdef int i
    cdef double result = 0.0

    for i in range(size):
        result += u[i] * v[i]

    return result


cdef double internal_vector_norm_(double *u, int norm, int size):
    """
    "Pure" C definition of vector norm for Cython compiler to be only used in the backend.
    
    Args:
        u: 
        norm: 
        size: 

    Returns:

    """

    cdef int i
    cdef double result = 0.0
    cdef double inverse_norm = 1.0 / norm

    for i in range(size):
        result += pow_C(u[i], norm)

    result = pow_C(result, inverse_norm)

    return result


cpdef double dot_product_(list u, list v):
    """
    Python API to use Cython dot_product implementation.
    
    Args:
        u list: Python list representing vector u  
        v list: Python list representing vector v

    Returns: Dot product of u and v

    """

    cdef int size
    if len(u) == len(v):
        size = len(u)
    else:
        raise ValueError("Vector size mismatch")

    cdef double *u_ = get_C_double_array_pointers(u, size)
    cdef double *v_ = get_C_double_array_pointers(v, size)

    cdef double result = internal_vector_dot_product_(u_, v_, size)

    release_C_pointer(u_)
    release_C_pointer(v_)

    return result
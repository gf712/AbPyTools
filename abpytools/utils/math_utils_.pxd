cdef double internal_vector_dot_product_(double *u, double *v, int size)
cdef double internal_vector_norm_(double *u, int norm, int size)

cdef class Matrix2D:
    cdef public int n_cols, n_rows
    cdef public list data
    cdef double* data_
    cdef list rows
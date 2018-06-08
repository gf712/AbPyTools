cdef class Matrix2D_backend:
    cdef public int n_cols, n_rows
    cdef values_
    cdef list data_
    cdef double* data_C
    cdef double** data_C_pointer
    cdef list rows

    cdef double* _get_row_pointer(self, int row)
    cdef double* _get_value_pointer(self, int row, int column)

cdef class Vector:
    cdef int size_
    cdef double* data_C

    @staticmethod
    cdef Vector create(double* ptr, int size)
    cdef double _get_element(self, int idx)
    cdef void _set_array_value(self, int idx, double value)
    cpdef double dot_product(self, Vector other)

# cdef class Vector:
#     cdef public list data
#     cdef double* data_
#     cdef public int size_
#
#     cpdef double dot_product(self, Vector other)

cdef double internal_vector_dot_product_(double *u, double *v, int size)
cdef double internal_vector_norm_(double *u, int norm, int size)

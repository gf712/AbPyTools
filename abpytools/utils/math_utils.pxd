cdef class Matrix2D_backend:
    cdef public int n_cols, n_rows, size_
    cdef list values_
    cdef double* data_C

    cdef double* _get_row_pointer(self, int row)
    cdef double* _get_value_pointer(self, int row, int column)
    cpdef list _get_row(self, int row)
    cpdef list _get_values(self)
    cdef double _get_value(self, int row, int col)
    cdef void _set_value(self, int row, int col, double value)

cdef class Vector:
    cdef int size_, derived_
    cdef double* data_C

    # additional instantiation methods
    @staticmethod
    cdef Vector create(double* ptr, int size)
    @staticmethod
    cdef Vector allocate(int size)
    @staticmethod
    cdef Vector allocate_with_value(int size, double value)

    cdef double _get_element(self, int idx)
    cdef void _set_array_value(self, int idx, double value)
    cpdef double dot_product(self, Vector other)
    cpdef Vector add(self, Vector other)
    cpdef Vector subtract(self, Vector other)
    cpdef Vector multiply(self, Vector other)
    cpdef Vector exp(self)
    cpdef double norm(self, int p=*)


cdef double internal_vector_dot_product_(double *u, double *v, int size)
cdef double internal_vector_norm_(double *u, int norm, int size)

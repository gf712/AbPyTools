import cython

ctypedef fused scalar_or_char:
    cython.floating
    cython.integral
    cython.char

cdef double* get_C_double_array_pointers(list a, int size)
cdef char* get_C_char_array_pointers(str a, int size)
cdef double* get_array_from_ptr(double* ptr, int size)
cdef void release_C_pointer(scalar_or_char *a)
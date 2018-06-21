import cython

cdef extern from "<malloc.h>" nogil:
    void *memalign(size_t alignment, size_t size);

ctypedef fused scalar_or_char:
    cython.floating
    cython.integral
    cython.char

ctypedef fused scalar:
    cython.floating
    cython.integral

cdef double* get_C_double_array_pointers(list a, int size)
cdef double** get_C_double_array_pp(double* a, int size)
cdef char* get_C_char_array_pointers(str a, int size)
cdef double** get_pp_from_ptr(double* ptr, int size)
cdef void get_array_from_ptr(double* ptr, double* a_, int size)
cdef void release_C_pointer(scalar_or_char *a)
cdef void release_C_pp(scalar_or_char** a)
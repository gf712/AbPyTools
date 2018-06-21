from libc.stdlib cimport malloc, free
from libc.string cimport strcpy
from cython.operator cimport dereference, postincrement
# from libc.stdio cimport printf


cdef double* get_C_double_array_pointers(list a, int size):

    """
    Internal function to convert a python list to C array
    :param a: 
    :param size: 
    :return: 
    """

    # allocate memory of a_ and b_ C arrays that will contain a copy of a and b lists, respectively
    IF SSE4_2:
        cdef double *a_ = <double *> memalign(16, size*sizeof(double *))
    ELSE:
        cdef double *a_ = <double *> malloc(size*sizeof(double *))


    if a_ is NULL:
        raise MemoryError("Failed to allocate memory!")

    cdef int i
    for i in range(size):
        a_[i] = a[i]
        i+=1

    return a_


cdef double** get_C_double_array_pp(double* a, int size):

    """
    Internal function to convert a python list to a C array of pointers that point to pointers (double**)
    :param a: 
    :param size: 
    :return: 
    """

    # allocate memory of a_ C arrays that will contain a copy of list a
    IF SSE4_2:
        cdef double **a_ = <double **> memalign(16, size*sizeof(double*))
    ELSE:
        cdef double **a_ = <double **> malloc(16, size*sizeof(double*))

    if a_ is NULL:
        raise MemoryError("Failed to allocate memory!")

    for i in range(size):
        a_[i] = a
        # printf("%f, %p\n", dereference(a_[i]), a_[i]);
        # print(i, ref, dereference(a_[i]), hex(<unsigned long>&a_[i]))
        postincrement(a)

    return a_


cdef char* get_C_char_array_pointers(str a, int size):

    """
    Internal function to convert a python string to C array
    :param a: 
    :param size: 
    :return: 
    """

    # allocate memory of a_ and b_ C arrays that will contain a copy of a and b lists, respectively
    IF SSE4_2:
        cdef char *a_ = <char*> memalign(16, size*sizeof(char *))
    ELSE:
        cdef char *a_ = <char*> malloc(16, size*sizeof(char *))

    if a_ is NULL:
        raise MemoryError("Failed to allocate memory!")

    cdef bytes a_py_bytes = a.encode('UTF-8')

    strcpy(a_, a_py_bytes)

    return a_


cdef void release_C_pointer(scalar_or_char *a):

    """
    Internal function to release memory allocated for pointer *a.
    :param a: 
    :return: 
    """

    free(a)


cdef void release_C_pp(scalar_or_char **a):

    free(a[0])
    free(a)


cdef double** get_pp_from_ptr(double* ptr, int size):
    """
    Internal function to get C array from a single ptr
    Args:
        ptr: C pointer 
        size: size of array

    Returns:

    """
    IF SSE4_2:
        cdef double **a_ = <double **> memalign(16, size*sizeof(double*))
    ELSE:
        cdef double **a_ = <double **> malloc(16, size*sizeof(double*))

    if a_ is NULL:
        raise MemoryError("Failed to allocate memory!")

    cdef int i
    for i in range(size):
        a_[i] = ptr
        # printf("%f, %d\n", dereference(a_[i]), a_[i])
        # print(i, ref, dereference(a_[i]), hex(<unsigned long>&a_[i]))
        postincrement(ptr)

    return a_


cdef void get_array_from_ptr(double* ptr, double* a_, int size):

    cdef int i

    for i in range(size):
        a_[0] = dereference(ptr)
        postincrement(ptr)
        postincrement(a_)

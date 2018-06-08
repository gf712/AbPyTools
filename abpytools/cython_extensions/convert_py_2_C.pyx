from libc.stdlib cimport malloc, free
from libc.string cimport strcpy
from cython.operator cimport dereference, postincrement


cdef double* get_C_double_array_pointers(list a, int size):

    """
    Internal function to convert a python list to C array
    :param a: 
    :param size: 
    :return: 
    """

    # allocate memory of a_ and b_ C arrays that will contain a copy of a and b lists, respectively
    cdef double *a_ = <double *> malloc(size*sizeof(double *))

    if a_ is NULL:
        raise MemoryError("Failed to allocate memory!")

    cdef int i
    for i in range(size):
        a_[i] = a[i]
        i+=1

    return a_


cdef char* get_C_char_array_pointers(str a, int size):

    """
    Internal function to convert a python string to C array
    :param a: 
    :param size: 
    :return: 
    """

    # allocate memory of a_ and b_ C arrays that will contain a copy of a and b lists, respectively
    cdef char *a_ = <char*> malloc(size*sizeof(char *))

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


cdef double* get_array_from_ptr(double* ptr, int size):
    """
    Internal function to get C array from a single ptr
    Args:
        ptr: C pointer 
        size: size of array

    Returns:

    """
    cdef double *a_ = <double *> malloc(size*sizeof(double *))

    if a_ is NULL:
        raise MemoryError("Failed to allocate memory!")

    cdef int i
    for i in range(size):
        a_[i] = dereference(ptr)
        postincrement(ptr)

    return a_

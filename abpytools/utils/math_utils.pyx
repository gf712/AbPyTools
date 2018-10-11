from libc.math cimport pow as pow_C
from abpytools.cython_extensions.convert_py_2_C cimport (get_C_double_array_pointers, release_C_pointer, memalign, posix_memalign)
import itertools
from cython.operator cimport dereference, postincrement
from libc.stdlib cimport malloc, free

cdef extern from "ops.h":
    double norm_op(double *A, int p, int size)
    void exp_op(double *A, double *B, int size)
    double l2_norm_op(double *A, int size)
    double l1_norm_op(const double *A, int size)
    # elementwise ops
    void subtract_op(double *A, double *B, double *C, int size)
    void multiply_op(double *A, double *B, double *C, int size)
    void add_op(double *A, double *B, double *C, int size)


cdef class NumericalBaseStructure:
    pass


cdef class Matrix2D_backend:
    """
    Cython class to manipulate double precision 2D matrix
    """

    def __init__(self, values_, ptrs=None, shape=None):
        # if ptr and shape are defined, values will be ignored and the data will be read from ptrs
        # check all rows have the same size
        self.values_ = values_
        self.n_rows = len(self.values_)
        self.n_cols = len(self.values_[0])
        self.size_ = self.n_cols * self.n_rows

        cdef int i
        for i in range(1, self.n_rows):
           if len(self.values_[i]) != self.n_cols:
               raise ValueError("All rows must have the same number of elements")

        # copy data to a contiguous C array
        self.data_C = get_C_double_array_pointers(list(itertools.chain(*self.values_)), self.size_)

    def _check_index(self, int row, int column):
        if row * self.n_cols + column < self.size_:
            raise IndexError("Requested element is outside array bounds")

    cdef double* _get_value_pointer(self, int row, int col):
        cdef double* ptr = &self.data_C[self.n_cols * row + col]

    cdef double* _get_row_pointer(self, int row):
        """
        Get pointer of first element of this row
        Args:
            row:

        Returns:

        """
        # pointer to first element of C array row
        cdef double* ptr = &self.data_C[self.n_cols * row]
        # print(dereference(ptr), hex(<unsigned long>&ptr))
        return ptr

    cpdef list _get_row(self, int row):
        result = []
        cdef int i
        cdef double* ptr = self._get_row_pointer(row)
        for i in range(self.n_cols):
            result.append(dereference(ptr))
            postincrement(ptr)

        return result

    def __dealloc__(self):
        # releases C memory allocation
        # print("Releasing Matrix memory ({})".format(id(self)))
        release_C_pointer(self.data_C)
        # release_C_pp(self.data_C)

    cdef double _get_value(self, int row, int col):
        return dereference(self._get_row_pointer(row)+col)

    cdef void _set_value(self, int row, int col, double value):
        cdef double* ptr =  self._get_row_pointer(row) + col
        ptr[0] = value


    def __getitem__(self, idx):
        if isinstance(idx, tuple) and len(idx) == 2:
            # return a single value (only supports rows and columns)
            # NOTE that there are no bound checks, so this can lead to odd behaviour,
            # e.g. asking for a[0, a.shape[1]+1] will return the first element of the second row
            # however this avoids wasting resources creating a new Vector object to get a single element
            row, column = idx
            return self._get_value(row, column)
        else:
            # create vector using data_ pointers
            return Vector.create(self._get_row_pointer(idx), self.n_cols)

    def __setitem__(self, idx, value):
        if isinstance(idx, tuple) and len(idx) == 2:
            row, column = idx
            self._set_value(row, column, value)
        else:
            raise NotImplementedError("Cannot assign vectors to rows in Matrix2D yet!")


    cpdef list _get_values(self):
        """
        Returns all values of C array
        Returns:

        """
        cdef list result = []
        cdef int i
        cdef int j
        cdef double* ptr
        cdef list row_result

        for i in range(self.n_rows):
            row_result = []
            ptr = self._get_row_pointer(i)  # this is more of a sanity check than anything else
            for j in range(self.n_cols):
                row_result.append(dereference(ptr))  # same as *ptr in C
                postincrement(ptr) # same as ptr++ in C
            result.append(row_result)

        return result



class Matrix2D(Matrix2D_backend):

    def __init__(self, values):
        """
        Lightweight numpy-like class for fast numerical calculations with Cython
        Args:
            values:
        """
        super().__init__(values)

    @property
    def values(self):
        """
        Returns data as a list of lists
        """
        return self._get_values()

    @property
    def shape(self):
        return self.n_rows, self.n_cols

    def __str__(self):

        final_string = "[{}]"
        intermediate_string = ""
        for i in range(self.n_rows):
            intermediate_string += "[{}]".format(', '.join([str(x) for x in self._get_row(i)]))

        return final_string.format(intermediate_string)

    def __repr__(self):

        repr_string = "array"
        first=True
        for x in self.__str__().split('\n'):
            if first:
                repr_string += "{}\n".format(x)
                first=False
            else:
                repr_string+="     {}\n".format(x)
        return repr_string


cdef class Vector:

    """
    Lightweight Cython vector class implementation that stores data as C arrays and uses these to perform calculations
    """

    def __init__(self, values_=None):
        if values_ is not None:
            self.size_ = len(values_)
            self.data_C = get_C_double_array_pointers(values_, self.size_)
            self.derived_=0


    @staticmethod
    cdef Vector create(double* ptr, int size):
        """
        Create a Vector from a pointer.
        
        Args:
            ptr: 
            size: 

        Returns:

        """
        cdef Vector vec = Vector()
        vec.size_ = size
        vec.data_C = ptr
        vec.derived_=1
        return vec

    @staticmethod
    cdef Vector allocate(int size):
        """
        
        Args:  
            size: 

        Returns:

        """
        IF SSE4_2:
            cdef double *value_pointers_ = <double *> memalign(16, size*sizeof(double))
        ELSE:
            cdef double *value_pointers_ = <double *> malloc(size*sizeof(double))

        if value_pointers_ is NULL:
            raise MemoryError("Failed to allocated memory")

        cdef Vector vec = Vector()
        vec.size_=size
        vec.data_C = value_pointers_
        vec.derived_=0

        return vec

    @staticmethod
    cdef Vector allocate_with_value(int size, double value):
        cdef Vector vec = Vector.allocate(size)
        cdef double *value_ptr = vec.data_C
        for i in range(size):
            value_ptr[0] = value
            postincrement(value_ptr)
        vec.derived_=0
        return vec


    cpdef double dot_product(self, Vector other):

        """
        Vector dot product.
        Args:
            other (Vector): a Vector object to calculate the dot product with

        Returns: the dot product between two vectors

        """

        if self.size_ != other.size:
            raise ValueError("Vector size mismatch")

        return internal_vector_dot_product_(self.data_C, other.data_C, self.size_)

    cdef double _get_element(self, int idx):
        return self.data_C[idx]

    cdef void _set_array_value(self, int idx, double value):
        cdef double* ptr = &self.data_C[idx]
        ptr[0] = value

    cpdef Vector exp(self):
        cdef Vector vec = Vector.allocate(self.size_)
        exp_op(self.data_C, vec.data_C, self.size_)
        return vec

    cpdef double norm(self, int p=2):
        if p==2:
            return l2_norm_op(self.data_C, self.size_)
        elif p == 1:
            return l1_norm_op(self.data_C, self.size_)
        else:
            return internal_vector_norm_(self.data_C, p, self.size_)

    @property
    def values(self):
        """
        Returns:
        """
        return [self._get_element(i) for i in range(self.size_)]

    def __dealloc__(self):
        if not self.derived_:
            if self.data_C != NULL:
                free(self.data_C)

    @property
    def size(self):
        return self.size_

    def __getitem__(self, int idx):
        return self._get_element(idx)

    def __setitem__(self, idx, value):
        self._set_array_value(idx, value)

    cpdef Vector subtract(self, Vector other):
        """
        Elementwise subtraction between two vectors.
        
        Args:
            other: 

        Returns:

        """

        # create empty vector
        cdef Vector vec = Vector.allocate(self.size_)

        subtract_op(self.data_C, other.data_C, vec.data_C, self.size_)

        return vec


    cpdef Vector add(self, Vector other):
        """
        Elementwise addition between two vectors.
        
        Args:
            other: 
    
        Returns:
    
        """

        # create empty vector
        cdef Vector vec = Vector.allocate(self.size_)

        add_op(self.data_C, other.data_C, vec.data_C, self.size_)

        return vec


    cpdef Vector multiply(self, Vector other):
        """
        Elementwise multiplication between two vectors.
        
        Args:
            other: 

        Returns:

        """
        # create empty vector
        cdef Vector vec = Vector.allocate(self.size_)

        multiply_op(self.data_C, other.data_C, vec.data_C, self.size_)

        return vec

    def __add__(self, other):
        return self.add(other)

    def __sub__(self, other):
        return self.subtract(other)

    def __mul__(self, other):
        return self.multiply(other)

    def __iter__(self):
        for x in range(self.size_):
            yield self[x]


cdef double internal_vector_dot_product_(double *u, double *v, int size):
    """
    "Pure" C definition of dot product for Cython compiler (to be only used in the backend).
    
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
    "Pure" C definition of vector norm for Cython compiler (to be only used in the backend).
    
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

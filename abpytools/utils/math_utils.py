def vector(A, B):
    """
    Creates a vector formed by points A and B
    :param A:
    :param B:
    :return:
    """

    try:
        _check_arrays(arrays=[A, B], check='type')
    except TypeError:
        if not isinstance(A, list) or not isinstance(B, list):
            raise TypeError("Expected A and B to be lists, but got {} and {} instead".format(type(A), type(B)))

    try:
        _check_arrays(arrays=[A, B], check='dim_min_size', min_size=1)
    except ValueError:
        raise ValueError("Both points must have at least 2 dimensions, "
                         "instead got list with {} and {} dimensions".format(len(A), len(B)))

    try:
        _check_arrays(arrays=[A, B], check='dims_consistent')
    except ValueError:
        raise ValueError("Expected A and B to have the same dimensions, "
                         "instead got dimensions {} and {}".format(len(A), len(B)))

    return [a - b for a, b in zip(A, B)]


def dot_product(u, v):
    """
    Calculates the dot product of two vectors u and v
    :param u:
    :param v:
    :return:
    """

    try:
        _check_arrays(arrays=[u, v], check='type')
    except TypeError:
        if not isinstance(u, list) or not isinstance(v, list):
            raise TypeError("Expected u and v to be lists, but got {} and {} instead".format(type(u), type(v)))
    try:
        _check_arrays(arrays=[u, v], check='dim_min_size', min_size=1)
    except ValueError:
        raise ValueError("Both vectors must have at least 2 dimensions, "
                         "instead got list with {} and {} dimensions".format(len(u), len(v)))

    try:
        _check_arrays(arrays=[u, v], check='dims_consistent')
    except ValueError:
        raise ValueError("Expected u and v to have the same dimensions, "
                         "instead got dimensions {} and {}".format(len(u), len(v)))

    return sum([u_i * v_i for u_i, v_i in zip(u, v)])


def magnitude(u):
    """
    Calculates the magnitude (euclidean norm) of vector u
    :param u: 
    :return: 
    """
    return sum([x ** 2 for x in u]) ** 0.5


def _check_arrays(arrays, check='dims_consistent', min_size=0):

    if check == 'dims_consistent':
        array_lengths = [len(x) for x in arrays]
        if len(set(array_lengths)) > 1:
            raise ValueError("Check array dimensions.")

    elif check == 'type':
        array_types = [type(x) for x in arrays]
        if len(set(array_types)) > 1 and list not in array_types:
            raise TypeError("Check the array list types.")

    elif check == 'dim_min_size':
        array_lengths = [len(x) for x in arrays]
        if len(set(array_lengths)) > 1:
            raise ValueError("Check array dimensions.")
        if list(set(array_lengths))[0] < min_size:
            raise ValueError("Not all elements are of size {}".format(min_size))

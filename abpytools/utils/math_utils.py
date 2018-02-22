def vector(A, B):
    """
    Creates a vector formed by points A and B
    :param A:
    :param B:
    :return:
    """

    _check_arrays(arrays=[A, B], checks=['type', 'dim_min_size', 'dims_consistent'],
                  min_size=1)

    return [a - b for a, b in zip(A, B)]


def dot_product(u, v):
    """
    Calculates the dot product of two vectors u and v
    :param u:
    :param v:
    :return:
    """

    _check_arrays(arrays=[u, v], checks=['type', 'dim_min_size', 'dims_consistent'],
                  min_size=1)

    return sum([u_i * v_i for u_i, v_i in zip(u, v)])


def magnitude(u):
    """
    Calculates the magnitude (euclidean norm) of vector u
    :param u: list representing a vector
    :return: float magnitude of vector
    """
    return sum([x ** 2 for x in u]) ** 0.5


def _check_arrays(arrays, checks=None, min_size=0):

    if checks is None:
        checks = ['dims_consistent']

    for check in checks:
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

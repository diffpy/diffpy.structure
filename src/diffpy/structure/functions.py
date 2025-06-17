import numpy as np


def dot_product(a, b):
    """Compute the dot product of two vectors of any size.

    Ensure that the inputs, a and b, are of the same size.
    The supported types are "array_like" objects, which can
    be converted to a NumPy array. Examples include lists and tuples.

    Parameters
    ----------
    a : array_like
        The first input vector.
    b : array_like
        The second input vector.

    Returns
    -------
    float
        The dot product of the two vectors.

    Examples
    --------
    Compute the dot product of two lists:
    >>> a = [1, 2, 3]
    >>> b = [4, 5, 6]
    >>> dot_product(a, b)
    32.0
    """
    return float(np.dot(a, b))

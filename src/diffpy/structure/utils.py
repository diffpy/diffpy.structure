#!/usr/bin/env python
##############################################################################
#
# diffpy.structure  by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2006 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################

"""Small shared functions.
"""

from collections.abc import Iterable as _Iterable

import numpy


def isiterable(obj):
    """``True`` if argument is iterable."""
    rv = isinstance(obj, _Iterable)
    return rv


def isfloat(s):
    """``True`` if argument can be converted to float."""
    try:
        float(s)
        return True
    except ValueError:
        pass
    return False


def atomBareSymbol(smbl):
    """Remove atom type string stripped of isotope and ion charge symbols.

    This function removes any blank, isotope numbers (0-9), leading hyphens (-), and ion charge
    symbols (1-9)(+-) from the given atom type string, returning only the bare element symbol.

    Parameters
    ----------
    smbl : str
        Atom type string that may include isotope numbers, ion charges, or hyphens.

    Returns
    -------
    str
        The bare element symbol.

    Examples
    --------
    >>> atomBareSymbol("Cl-")
    'Cl'
    >>> atomBareSymbol("Ca2+")
    'Ca'
    >>> atomBareSymbol("12-C")
    'C'
    """

    rv = smbl.strip().lstrip("0123456789-").rstrip("123456789+-")
    return rv


# Helpers for the Structure class --------------------------------------------


def _linkAtomAttribute(attrname, doc, toarray=numpy.array):
    """Create property wrapper that maps the specified atom attribute.

    The returned property object provides convenient access to atom
    attributes from the owner `Structure` class.

    Parameters
    ----------
    attrname : str
        The string name of the `Atom` class attribute to be mapped.
    doc : str
        The docstring for the property wrapper.
    toarray : callable, Optional
        Factory function that converts list of attributes to `numpy.ndarray`.
        Use `numpy.char.array` for string attributes.

    Return a property object.
    """
    from itertools import repeat
    from operator import setitem

    _all = slice(None)

    def fget(self):
        va = toarray([getattr(a, attrname) for a in self])
        return va

    def fset(self, value):
        n = len(self)
        if n == 0:
            return
        v0 = getattr(self[0], attrname)
        # replace scalar values, but change array attributes in place
        if numpy.isscalar(v0):

            def setvalue(a, v):
                return setattr(a, attrname, v)

        else:

            def setvalue(a, v):
                return setitem(getattr(a, attrname), _all, v)

        # avoid broadcasting if the new value is a scalar
        if numpy.isscalar(value):
            genvalues = repeat(value)
        else:
            genvalues = numpy.broadcast_to(value, (n,) + numpy.shape(v0))
        for a, v in zip(self, genvalues):
            setvalue(a, v)
        return

    rv = property(fget, fset, doc=doc)
    return rv

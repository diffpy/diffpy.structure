#!/usr/bin/env python
##############################################################################
#
# diffpy.Structure  by DANSE Diffraction group
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

import numpy


def isfloat(s):
    """True if argument can be converted to float"""
    try:
        float(s)
        return True
    except ValueError:
        pass
    return False


def atomBareSymbol(smbl):
    '''Remove atom type string stripped of isotope and ion charge symbols.
    This removes blank and leading [0-9-] or trailing [1-9][+-] characters.

    smbl -- atom type string such as "Cl-", "Ca2+" or "12-C".

    Return bare element symbol.
    '''
    rv = smbl.strip().lstrip('0123456789-').rstrip('123456789+-')
    return rv

# Helpers for the Structure class --------------------------------------------

def _linkAtomAttribute(attrname, doc, toarray=numpy.array):
    '''Create property wrapper that maps the specified atom attribute.
    The returned property object provides convenient access to atom
    attributes from the owner Structure class.

    attrname -- string name of the Atom class attribute to be mapped
    doc      -- docstring of the property wrapper
    toarray  -- factory function that converts list of attributes to
                numpy.array.  Use numpy.char.array for string attributes.

    Return a property object.
    '''
    def fget(self):
        va = toarray([getattr(a, attrname) for a in self])
        return va
    def fset(self, value):
        if len(self) == 0:  return
        # dummy array va helps to broadcast the value to proper iterable
        va = numpy.asarray(len(self) * [getattr(self[0], attrname)])
        for a, v in zip(self, numpy.broadcast_arrays(va, value)[1]):
            setattr(a, attrname, v)
        return
    rv = property(fget, fset, doc=doc)
    return rv

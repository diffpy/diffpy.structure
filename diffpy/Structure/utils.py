##############################################################################
#
# Structure         by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2006 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################

"""Small shared functions.
"""

__id__ = "$Id$"

import types
import numpy


def isfloat(s):
    """True if argument can be converted to float"""
    try:
        x = float(s)
        return True
    except ValueError:
        pass
    return False

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
    rv = property(fget, fset, doc)
    return rv

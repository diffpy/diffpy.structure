#!/usr/bin/env python
##############################################################################
#
# diffpy.structure  by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2008 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Chris Farrow, Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################

"""This module contains functions for simple `Structure` manipulation.
"""

import numpy

from diffpy.structure import Atom, Structure


def supercell(S, mno):
    """
    Perform supercell expansion for a `Structure`.

    New `lattice` parameters are multiplied and fractional coordinates
    divided by corresponding multiplier. New `Atoms` are grouped with
    their source in the original cell.

    Parameters
    ----------
    S : Structure
        An instance of `Structure` from `diffpy.structure`.
    mno : array_like
        Sequence of 3 integers for cell multipliers along
        the `a`, `b` and `c` axes.

    Returns
    -------
    Structure
        A new `Structure` instance representing the expanded supercell.

    Raises
    ------
    TypeError
        `S` is not a `Structure` instance.
    ValueError
        Invalid `mno` argument.
    """
    # check arguments
    if len(mno) != 3:
        emsg = "Argument mno must contain 3 numbers."
        raise ValueError(emsg)
    elif min(mno) < 1:
        emsg = "Multipliers must be greater or equal 1"
        raise ValueError(emsg)
    if not isinstance(S, Structure):
        emsg = "The first argument must be a Structure instance."
        raise TypeError(emsg)

    # convert mno to a tuple of integers so it can be used as range limit.
    mno = (int(mno[0]), int(mno[1]), int(mno[2]))

    # create return instance
    newS = Structure(S)
    if mno == (1, 1, 1):
        return newS

    # back to business
    ijklist = [(i, j, k) for i in range(mno[0]) for j in range(mno[1]) for k in range(mno[2])]
    # numpy.floor returns float array
    mnofloats = numpy.array(mno, dtype=float)

    # build a list of new atoms
    newAtoms = []
    for a in S:
        for ijk in ijklist:
            adup = Atom(a)
            adup.xyz = (a.xyz + ijk) / mnofloats
            newAtoms.append(adup)
    # newS can own references in newAtoms, no need to make copies
    newS.__setitem__(slice(None), newAtoms, copy=False)

    # take care of lattice parameters
    newS.lattice.setLatPar(a=mno[0] * S.lattice.a, b=mno[1] * S.lattice.b, c=mno[2] * S.lattice.c)
    return newS


# End of supercell

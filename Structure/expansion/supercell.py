#!/usr/bin/env python
"""This module contains methods for simple structure manipulation."""

__id__ = "$Id$"

from diffpy.Structure import Structure, Atom
import numpy


def createSuperCell(S, l = 1, m = 1, n = 1):
    """Perform supercell expansion for a structure.
    
    New lattice parameters are multiplied and fractional coordinates divided by
    corresponding multiplier.  New atoms are grouped with their source in the
    original cell.

    Arguments
    S       --  A Structure instance
    l       --  Cell multiplier along a-axis (positive integer, default 1)
    m       --  Cell multiplier along b-axis (positive integer, default 1)
    n       --  Cell multiplier along c-axis (positive integer, default 1)
    
    Returns a new expanded structure instance
    """
    # check arguments
    lmn = tuple(map(int, (l,m,n)))
    if min(lmn) < 1:
        raise ValueError("Multipliers must be >= 1")
    if not isinstance(S, Structure):
        raise TypeError("Must pass a structure instance as the first argument")

    newS = Structure(S)
    if lmn == (1,1,1):
        return newS

    # back to business
    ijklist = [(i,j,k) 
                for i in range(lmn[0]) 
                    for j in range(lmn[1]) 
                        for k in range(lmn[2])]
    lmnfloats = numpy.array(lmn[:], dtype=float)

    # build a list of new atoms
    newAtoms = []
    for a in S:
        for ijk in ijklist:
            adup = Atom(a)
            adup.xyz = (a.xyz + ijk)/lmnfloats
            newAtoms.append(adup)
    newS[:] = newAtoms

    # take care of lattice parameters
    newS.lattice.setLatPar(
            a=lmn[0]*S.lattice.a,
            b=lmn[1]*S.lattice.b,
            c=lmn[2]*S.lattice.c )
    return newS

if __name__ == "__main__":

    import os.path
    datadir = "../../tests/testdata"
    S = Structure()
    S.read(os.path.join(datadir, "Ni.stru"), "pdffit")
    newS = createSuperCell(S, 2, 2, 2)
    newS.write("Ni_2x2x2.stru", "pdffit")

    S = Structure()
    S.read(os.path.join(datadir, "CdSe-wurtzite.stru"), "pdffit")
    newS = createSuperCell(S, 2, 2, 2)
    newS.write("CdSe_2x2x2.stru", "pdffit")


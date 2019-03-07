#!/usr/bin/env python
##############################################################################
#
# diffpy.Structure  by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2008 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################

"""Utilities for making shapes."""

def findCenter(S):
    """Find the approximate center atom of a structure.

    The center of the structure is the atom closest to (0.5, 0.5, 0.5)

    Returns the index of the atom.
    """
    best = -1
    bestd = len(S)
    center = [0.5, 0.5, 0.5] # the cannonical center

    for i in range(len(S)):
        d = S.lattice.dist(S[i].xyz, center)
        if d < bestd:
            bestd = d
            best = i

    return best


if __name__ == "__main__":
    # FIXME ... remove or convert to unit test
    import os.path
    datadir = "../../tests/testdata"
    S = Structure()
    S.read(os.path.join(datadir, "CdSe_bulk.stru"), "pdffit")
    newS = makeEllipsoid(S, 20)
    newS.write("CdSe_d20.stru", "pdffit")
    newS = makeEllipsoid(S, 20, 10, 10)
    newS.write("CdSe_a20_b10_c10.stru", "pdffit")
    newS = makeEllipsoid(S, 20, 15, 10)
    newS.write("CdSe_a20_b15_c10.stru", "pdffit")
    S = Structure()
    S.read(os.path.join(datadir, "Ni.stru"), "pdffit")
    newS = makeEllipsoid(S, 20)
    newS.write("Ni_d20.stru", "pdffit")
    newS = makeEllipsoid(S, 20, 4)
    newS.write("Ni_a20_b4_c20.stru", "pdffit")
    newS = makeEllipsoid(S, 20, 15, 10)
    newS.write("Ni_a20_b15_c10.stru", "pdffit")

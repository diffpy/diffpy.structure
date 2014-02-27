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

"""Make a spheroid nanoparticle from a template structure."""

from math import ceil
from diffpy.Structure import Structure, Atom
from diffpy.Structure.expansion.shapeUtils import findCenter

def makeCuboctahedron(S, dist):
    """Make a cuboctahedron nanoparticle.

    Arguments
    S       --  A Structure instance
    dist    --  Distance from center to nearest face

    Returns a new structure instance
    """

    # Create a supercell large enough for the ellipsoid
    frac = S.lattice.fractional((dist, dist, dist))
    mno = map(ceil, 2*frac)
    # Make the supercell
    from diffpy.Structure.expansion import supercell
    newS = supercell(S, mno)
    lat = newS.lattice

    # Find the central atom
    ncenter = findCenter(newS)

    # Make the cuboctahedron template
    from geometry.composites import cuboctahedron
    from geometry.operations import translate, rotate
    c0 = translate(cuboctahedron(dist), lat.cartesian(newS[ncenter].xyz))

    # Cut out an octahedron
    from geometry import locate

    N = len(newS)
    j = N
    for i in xrange(N):

        xyz = lat.cartesian(newS[N-1-i].xyz)
        if locate(xyz, c0) == 1:
            newS.pop(N-1-i)

    return newS

if __name__ == "__main__":

    import os.path
    datadir = "../../../tests/testdata"
    S = Structure()
    S.read(os.path.join(datadir, "CdSe_bulk.stru"), "pdffit")
    newS = makeCuboctahedron(S, 12)
    newS.write("CdSe_cuboct24.stru", "pdffit")
    S = Structure()
    S.read(os.path.join(datadir, "Ni.stru"), "pdffit")
    newS = makeCuboctahedron(S, 10)
    newS.write("Ni_cuboct20.stru", "pdffit")

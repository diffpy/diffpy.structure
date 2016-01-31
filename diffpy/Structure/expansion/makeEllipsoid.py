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
from numpy import array
from diffpy.Structure import Structure
from diffpy.Structure.expansion.shapeUtils import findCenter

def makeSphere(S, radius):
    """Create a spherical nanoparticle.

    Arguments
    S       --  A Structure instance
    radius  --  primary equatorial radius (along x-axis)

    Returns a new structure instance
    """
    return makeEllipsoid(S, radius)

def makeEllipsoid(S, a, b=None, c=None):
    """Cut a structure out of another one.

    Arguments
    S       --  A Structure instance
    a       --  primary equatorial radius (along x-axis)
    b       --  secondary equatorial radius (along y-axis). If b is None
                (default) then it is set equal to a
    c       --  polar radius (along z-axis). If c is None (default), then it is
                set equal to a.

    Returns a new structure instance
    """
    if b is None: b = a
    if c is None: c = a
    sabc = array([a, b, c])

    # Create a supercell large enough for the ellipsoid
    frac = S.lattice.fractional(sabc)
    mno = map(ceil, 2*frac)
    mno = max(map(ceil, 2*frac))*array([1,1,1])
    # Make the supercell
    from diffpy.Structure.expansion import supercell
    newS = supercell(S, mno)
    lat = newS.lattice

    # Find the central atom
    ncenter = findCenter(newS)

    cxyz = lat.cartesian(newS[ncenter].xyz)

    delList = []
    N = len(newS)
    j = N
    for i in xrange(N):
        j -= 1

        # Calculate (x/a)**2 + (y/b)**2 + (z/c)**2
        xyz = lat.cartesian(newS[j].xyz)
        darray = ((xyz-cxyz)/sabc)**2
        d = sum(darray)**0.5

        # Discard atom if (x/a)**2 + (y/b)**2 + (z/c)**2 > 1
        if d > 1:
            delList.append(j)

    for i in delList:
        newS.pop(i)

    return newS


if __name__ == "__main__":
    import os.path
    datadir = "../../tests/testdata"
    S = Structure()
    S.read(os.path.join(datadir, "CdSe_bulk.stru"), "pdffit")
    newS = makeEllipsoid(S, 12)
    newS.write("CdSe_d24.stru", "pdffit")
    newS = makeEllipsoid(S, 20, 10, 10)
    newS.write("CdSe_a20_b10_c10.stru", "pdffit")
    newS = makeEllipsoid(S, 20, 15, 10)
    newS.write("CdSe_a20_b15_c10.stru", "pdffit")
    S = Structure()
    S.read(os.path.join(datadir, "Ni.stru"), "pdffit")
    newS = makeEllipsoid(S, 10)
    newS.write("Ni_d20.stru", "pdffit")
    newS = makeEllipsoid(S, 20, 4)
    newS.write("Ni_a20_b4_c20.stru", "pdffit")
    newS = makeEllipsoid(S, 20, 15, 10)
    newS.write("Ni_a20_b15_c10.stru", "pdffit")

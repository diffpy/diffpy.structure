#!/usr/bin/env python
##############################################################################
#
# diffpy.structure  by DANSE Diffraction group
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

from diffpy.structure import Structure
from diffpy.structure.expansion.shapeutils import findCenter


def makeSphere(S, radius):
    """Create a spherical nanoparticle.

    Parameters
    ----------
    S : Structure
        A `Structure` instance.
    radius : float
        Primary equatorial radius (along x-axis).

    Returns
    -------
    Structure
        A new `Structure` instance.
    """
    return makeEllipsoid(S, radius)


def makeEllipsoid(S, a, b=None, c=None):
    """
    Cut a `Structure` out of another one.

    Parameters
    ----------
    S : Structure
        A `Structure` instance.
    a : float
        Primary equatorial radius (along x-axis).
    b : float, Optional
        Secondary equatorial radius (along y-axis). If `b` is ``None``
        (default), then it is set equal to `a`.
    c : float, Optional
        Polar radius (along z-axis). If `c` is ``None`` (default), then it is
        set equal to `a`.

    Returns
    -------
    Structure :
        A new `Structure` instance.
    """
    if b is None:
        b = a
    if c is None:
        c = a
    sabc = array([a, b, c])

    # Create a supercell large enough for the ellipsoid
    frac = S.lattice.fractional(sabc)
    # FIXME - this looks fishy for non-orthogonal lattices
    mno = max(ceil(2 * xi) for xi in frac) * array([1, 1, 1])
    # Make the supercell
    from diffpy.structure.expansion import supercell

    newS = supercell(S, mno)
    lat = newS.lattice

    # Find the central atom
    ncenter = findCenter(newS)

    cxyz = lat.cartesian(newS[ncenter].xyz)

    delList = []
    N = len(newS)
    j = N
    for i in range(N):
        j -= 1

        # Calculate (x/a)**2 + (y/b)**2 + (z/c)**2
        xyz = lat.cartesian(newS[j].xyz)
        darray = ((xyz - cxyz) / sabc) ** 2
        d = sum(darray) ** 0.5

        # Discard atom if (x/a)**2 + (y/b)**2 + (z/c)**2 > 1
        if d > 1:
            delList.append(j)

    for i in delList:
        newS.pop(i)

    return newS


# ----------------------------------------------------------------------------

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

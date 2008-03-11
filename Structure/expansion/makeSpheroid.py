#!/usr/bin/env python
"""Make a spheroid nanoparticle from a template structure."""

from diffpy.Structure import Structure, Atom
from numpy import array

def makeSpheroid(S, radius, v = 1.0):
    """Cut a structure out of another one.

    Arguments
    S       --  A Structure instance
    radius  --  radius of spheroid
    v       --  ellipticity of spheroid (defined along c-axis, default 1.0)

    Returns a new structure instance
    """

    abc = array([S.lattice.a, S.lattice.b, S.lattice.c])
    sabc = array([radius, radius, v*radius])

    # Create a supercell large enough for the ellipsoid
    l = ceil(radius/S.lattice.a) + 1
    m = ceil(radius/S.lattice.b) + 1
    n = ceil(v*radius/S.lattice.c) + 1
    # Make these odd for simplicty
    if not l%2: l+=1
    if not m%2: m+=1
    if not n%2: n+=1

    # Find the central atom
    cl = l/2 + 1
    cm = m/2 + 1
    cn = n/2 + 1
    ci = 

    # This will become our [0,0,0]
    cxyz += shift

    for i in range(len(S), 0, -1):
        i -= 1
        S[i].xyz -= cxyz
        S2[i].xyz -= cxyz

        darray = (S[i].xyz*abc/sabc)**2
        d = sum(darray)**0.5

        if d > 1:
            del S2[i]

    return S2


def findCenter(S):
    """find the approximate center atom of a structure.

    Returns the xyz coordinates of the atom.
    """

    for i in range(len(S)):
        if (S[i].xyz < 0.6).all() and (S[i].xyz > 0.4).all():
            return array(S[i].xyz)

    raise RuntimeError("Can't find it!")


if __name__ == "__main__":


    S = cutSpheroid("Ni_10x10x10.stru", 10.0, 1.0, [0,0,0.05])
    S.write("Ni_d20b.xcfg", "xcfg")
    S.write("Ni_d20b.stru", "pdffit")

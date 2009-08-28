#!/usr/bin/env python

"""Demonstrate ways of building structure from Cartesian Coordinates.
"""

from diffpy.Structure import Structure
from diffpy.Structure import Lattice
from diffpy.Structure import Atom


def listCoordinates(stru):
    """List fractional and Cartesian coordinates in a structure.

    stru    -- instande of diffpy.Structure

    No return value.
    """
    for a in stru:
        print "%-2s  xyz=(%.4f, %.4f, %.4f)  xyz_cartn=(%.4f, %.4f, %.4f)" % (
                (a.element, ) + tuple(a.xyz) + tuple(a.xyz_cartn))
    return


# Structure instance has lattice attribute which describes coordinate
# system for fractional coordinates.  It defaults to Cartesian system.

stru1 = Structure(atoms=[
    Atom('C', [0, 0, 0]),
    Atom('C', [1, 1, 1]),
    Atom('C', [1, 2, 3]),
    Atom('C', [3, 2, 1]),
    ])

print "stru1 uses the default coordinate system"
print "stru1.lattice.abcABG()=" + str(stru1.lattice.abcABG())
listCoordinates(stru1)
print

# The placeInLattice method sets a new coordinate system while preserving
# the same Cartesian positions of all atoms.
stru2 = Structure(stru1)
lattice2 = Lattice(5, 6, 7, 60, 70, 80)
stru2.placeInLattice(lattice2)

print "stru2 is a copy of stru1 placed in differenc lattice"
print "stru2.lattice.abcABG()=" + str(stru2.lattice.abcABG())
listCoordinates(stru2)
print

# Finally to place atom at a given Cartesian position, one can
# set its xyz_cartn attribute.

lattice3 = Lattice(6, 8, 9, 90, 90, 90)
stru3 = Structure(lattice=lattice3)
# add 4 carbon atoms
stru3.addNewAtom("C")
stru3.getLastAtom().xyz_cartn = (0, 0, 0)
stru3.addNewAtom("C")
stru3[-1].xyz_cartn = (1, 1, 1)
stru3.addNewAtom("C")
stru3.getLastAtom().xyz_cartn = (1, 2, 3)
stru3.addNewAtom("C")
stru3.getLastAtom().xyz_cartn = (3, 2, 1)

print "stru3 atom coordinates were defined using xyz_cartn"
print "stru3.lattice.abcABG()=" + str(stru3.lattice.abcABG())
listCoordinates(stru3)
print

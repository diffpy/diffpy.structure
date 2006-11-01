########################################################################
#
# <PackageName>     by DANSE Diffraction group
#                   Simon J.L. Billinge
#                   Michigan State University
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See COPYRIGHT.txt for copying and usage conditions.
# See LICENSE.txt for license information.
#
########################################################################

"""class Atom for storing properties of a single atom"""

__id__ = "$Id$"

import numpy

class Atom:
    """Atom --> class for storing atom information

    Data members:
        element   -- type of the atom
        xyz       -- fractional coordinates
        occupancy -- fractional occupancy
        U         -- anisotropic thermal displacement coefficients
        name      -- atom label
    """

    def __init__(self, atype, xyz=None, name=None, occupancy=None, U=None):
        """Create atom of a specified type at given lattice coordinates.
        Atom(a) creates a copy of Atom instance a.

        atype       -- element symbol string or Atom instance
        xyz         -- fractional coordinates
        name        -- atom label
        occupancy   -- fractional occupancy
        U           -- anisotropic thermal displacement coefficients
        """
        # declare data members
        self.element = None
        self.xyz = numpy.zeros(3, dtype=float)
        self.name = ''
        self.occupancy = 1.0
        self.U = numpy.zeros((3,3), dtype=float)
        # assign them as needed
        if isinstance(atype, Atom):
            self.element = atype.element
            self.xyz = numpy.array(atype.xyz)
            self.name = atype.name
            self.occupancy = atype.occupancy
            self.U = numpy.array(atype.U)
        else:
            self.element = atype
        # take care of remaining arguments
        if xyz is not None:         self.xyz = numpy.array(xyz)
        if name is not None:        self.name = name
        if occupancy is not None:   self.occupancy = float(occupancy)
        if U is not None:           self.U = numpy.array(U)
        return

    def Uiso(self):
        """equivalent isotropic temperature factor"""
        Uiso = numpy.trace(self.U)/3.0
        return Uiso

    def Biso(self):
        """equivalent isotropic Debye-Waler temperature factor"""
        Biso = 8 * numpy.pi**2 * self.Uiso()
        return Biso

    def setUiso(self, Uiso):
        """set temperature coefficients matrix U to Uiso*identity(3)"""
        self.U = Uiso * numpy.identity(3, dtype=float)
        return

    def setBiso(self, Biso):
        """set matrix U to isotropic Debye-Waler factor Biso"""
        self.U = self.setUiso( Biso/(8.0*numpy.pi**2) )
        return

    def __repr__(self):
        """simple string representation"""
        s = "%-4s %8.6f %8.6f %8.6f %6.4f" % (self.element,
                self.xyz[0], self.xyz[1], self.xyz[2],
                self.occupancy)
        return s

# End of Atom

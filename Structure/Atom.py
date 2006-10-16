"""class Atom for storing properties of a single atom"""

__id__ = "$Id$"

import numpy as num

##############################################################################
class Atom:
    """Atom --> class for storing atom information

    Data members:
        element   -- type of the atom
        xyz       -- fractional coordinates
        occupancy -- fractional occupancy
        U         -- matrix of temperature coefficients
        name      -- atom label
    """

    def __init__(self, element, xyz, name='', occupancy=1.0, U=None):
        """create atom of a specified type at given lattice coordinates"""
        self.element = element
        self.xyz = num.array(xyz, dtype=float)
        self.name = name
        if U is None:
            self.U = num.zeros((3,3), dtype=float)
        else:
            self.U = num.array(U, dtype=float)
        self.occupancy = occupancy
        return

    def Uiso(self):
        """equivalent isotropic temperature factor"""
        Uiso = num.trace(self.U)/3.0
        return Uiso

    def Biso(self):
        """equivalent isotropic Debye-Waler temperature factor"""
        Biso = 8 * num.pi**2 * self.Uiso()
        return Biso

    def setUiso(self, Uiso):
        """set temperature coefficients matrix U to Uiso*identity(3)"""
        self.U = Uiso * num.identity(3, dtype=float)
        return

    def setBiso(self, Biso):
        """set matrix U to isotropic Debye-Waler factor Biso"""
        self.U = self.setUiso( Biso/(8.0*num.pi**2) )
        return

    def __repr__(self):
        """simple string representation"""
        s = "%-4s %8.6f %8.6f %8.6f %6.4f" % (self.element,
                self.xyz[0], self.xyz[1], self.xyz[2],
                self.occupancy)
        return s

# End of Atom

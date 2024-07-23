#!/usr/bin/env python
##############################################################################
#
# diffpy.structure  Complex Modeling Initiative
#                   (c) 2016 Brookhaven Science Associates,
#                   Brookhaven National Laboratory.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################

"""
Unit tests for the Atom class.
"""


import unittest

import numpy

from diffpy.structure.atom import Atom
from diffpy.structure.lattice import Lattice

# ----------------------------------------------------------------------------


class TestAtom(unittest.TestCase):

    def test___init__(self):
        """check Atom.__init__()"""
        a = Atom()
        self.assertEqual("", a.element)
        self.assertTrue((a.xyz == 0).all())
        self.assertEqual("", a.label)
        self.assertEqual(1.0, a.occupancy)
        self.assertFalse(a.anisotropy)
        self.assertTrue((a.U == 0).all())
        self.assertTrue(a.lattice is None)
        # check initialization with arguments
        a1 = Atom("C", xyz=(1, 2, 3), Uisoequiv=0.005)
        self.assertEqual("C", a1.element)
        self.assertTrue(numpy.array_equal([1, 2, 3], a1.xyz))
        self.assertFalse(a1.anisotropy)
        self.assertEqual(0.005, a1.Uisoequiv)
        # initialize with anisotropic displacement parameters
        uani = numpy.identity(3, dtype=float) * numpy.array([1, 2, 3]) * 0.01
        a2 = Atom("C", U=uani)
        self.assertTrue(numpy.array_equal(uani, a2.U))
        self.assertTrue(a2.anisotropy)
        a3 = Atom("C", Uisoequiv=0.02, anisotropy=True)
        self.assertTrue(a3.anisotropy)
        self.assertEqual(a3.U[2, 2], 0.02)
        self.assertRaises(ValueError, Atom, "C", Uisoequiv=0.02, U=uani)
        return

    #   def test_msdLat(self):
    #       """check Atom.msdLat()
    #       """
    #       return
    #
    #   def test_msdCart(self):
    #       """check Atom.msdCart()
    #       """
    #       return
    #
    #   def test___repr__(self):
    #       """check Atom.__repr__()
    #       """
    #       return
    #
    #   def test___copy__(self):
    #       """check Atom.__copy__()
    #       """
    #       return

    def test_xyz_cartn(self):
        """check Atom.xyz_cartn property"""
        hexagonal = Lattice(1, 1, 1, 90, 90, 120)
        a0 = Atom("C", [0, 0, 0], lattice=hexagonal)
        a1 = Atom("C", [1, 1, 1], lattice=hexagonal)
        self.assertTrue(all(a0.xyz_cartn == 0))
        rc1 = numpy.array([0.75**0.5, 0.5, 1])
        self.assertTrue(numpy.allclose(rc1, a1.xyz_cartn))
        a1.xyz_cartn[2] = 0
        self.assertTrue(numpy.allclose([1, 1, 0], a1.xyz))
        a1.xyz_cartn[:2] = 0
        self.assertTrue(all(a1.xyz == 0))
        a3 = Atom("C", [1, 2, 3])
        self.assertTrue(numpy.array_equal(a3.xyz, a3.xyz_cartn))
        a3.xyz_cartn = 1.3
        self.assertTrue(all(1.3 == a3.xyz_cartn))
        self.assertTrue(all(1.3 == a3.xyz))
        return


#   def test__get_anisotropy(self):
#       """check Atom._get_anisotropy()
#       """
#       return
#
#   def test__set_anisotropy(self):
#       """check Atom._set_anisotropy()
#       """
#       return
#
#   def test__get_U(self):
#       """check Atom._get_U()
#       """
#       return
#
#   def test__set_U(self):
#       """check Atom._set_U()
#       """
#       return
#
#   def test__get_Uij(self):
#       """check Atom._get_Uij()
#       """
#       return
#
#   def test__set_Uij(self):
#       """check Atom._set_Uij()
#       """
#       return
#
#   def test__get_Uisoequiv(self):
#       """check Atom._get_Uisoequiv()
#       """
#       return
#
#   def test__set_Uisoequiv(self):
#       """check Atom._set_Uisoequiv()
#       """
#       return
#
#   def test__get_Bisoequiv(self):
#       """check Atom._get_Bisoequiv()
#       """
#       return
#
#   def test__set_Bisoequiv(self):
#       """check Atom._set_Bisoequiv()
#       """
#       return

# End of class TestAtom

# ----------------------------------------------------------------------------

if __name__ == "__main__":
    unittest.main()

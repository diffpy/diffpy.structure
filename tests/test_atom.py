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
"""Unit tests for the Atom class."""


import unittest

import numpy
import pytest

from diffpy.structure.atom import Atom
from diffpy.structure.lattice import Lattice

# ----------------------------------------------------------------------------


class TestAtom(unittest.TestCase):

    def test___init__(self):
        """Check Atom.__init__()"""
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

    def test_msdLat(self):
        """Check Atom.msdLat."""
        hexagonal = Lattice(1, 1, 1, 90, 90, 120)
        atom_1 = Atom("C", [0, 0, 0], lattice=hexagonal, Uisoequiv=0.0123)
        assert atom_1.msdLat([1, 2, 3]) == pytest.approx(0.0123, rel=0, abs=1e-15)
        assert atom_1.msdLat([9, 0, -4]) == pytest.approx(0.0123, rel=0, abs=1e-15)
        U = numpy.array(
            [
                [0.010, 0.002, 0.001],
                [0.002, 0.020, 0.003],
                [0.001, 0.003, 0.030],
            ],
            dtype=float,
        )
        atom_2 = Atom("C", [0, 0, 0], lattice=hexagonal, U=U)

        vc = numpy.array([1.2, -0.3, 0.7], dtype=float)
        vl = hexagonal.fractional(vc)

        assert atom_2.msdLat(vl) == pytest.approx(atom_2.msd_cart(vc), rel=1e-13, abs=1e-13)

    def test_msdCart(self):
        """Check Atom.msdCart."""
        hexagonal = Lattice(1, 1, 1, 90, 90, 120)
        atom_1 = Atom("C", [0, 0, 0], lattice=hexagonal, Uisoequiv=0.0456)
        assert atom_1.msdCart([1, 0, 0]) == pytest.approx(0.0456, rel=0, abs=1e-15)
        assert atom_1.msdCart([0, 5, -2]) == pytest.approx(0.0456, rel=0, abs=1e-15)
        assert atom_1.msdCart([0, 5, -2]) == pytest.approx(0.0456, rel=0, abs=1e-15)

        U = numpy.array(
            [
                [0.011, 0.001, 0.000],
                [0.001, 0.019, 0.002],
                [0.000, 0.002, 0.027],
            ],
            dtype=float,
        )
        atom_2 = Atom("C", [0, 0, 0], lattice=hexagonal, U=U)

        vc = numpy.array([0.4, 1.1, -0.6], dtype=float)
        assert atom_2.msdCart(vc) == pytest.approx(atom_2.msdCart(3.7 * vc), rel=1e-13, abs=1e-13)

    def test_xyz_cartn(self):
        """Check Atom.xyz_cartn property."""
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
@pytest.mark.parametrize(
    "uiso, lattice_vector",
    [  # C1: isotropic displacement, msd is direction-independent in lattice coordinates.
        # Expected the msd_latt equals Uisoequiv for any direction.
        (0.0123, [1, 2, 3]),
    ],
)
def test_msd_latt_isotropic(uiso, lattice_vector):
    """Check Atom.msd_latt()."""
    hexagonal = Lattice(1, 1, 1, 90, 90, 120)
    atom = Atom("C", [0, 0, 0], lattice=hexagonal, Uisoequiv=uiso)
    actual = atom.msd_latt(lattice_vector)
    expected = pytest.approx(uiso, rel=0, abs=1e-15)
    assert actual == expected


@pytest.mark.parametrize(
    "U, cartesian_vector",
    [  # C2: anisotropic displacement with same physical direction expressed in lattice vs cartesian coords
        # Expected msd_latt(fractional(cartesian_vector)) == msd_cart(cartesian_vector)
        (
            numpy.array(
                [
                    [0.010, 0.002, 0.001],
                    [0.002, 0.020, 0.003],
                    [0.001, 0.003, 0.030],
                ],
                dtype=float,
            ),
            numpy.array([1.2, -0.3, 0.7], dtype=float),
        ),
        (
            numpy.array(
                [
                    [0.018, -0.001, 0.002],
                    [-0.001, 0.012, 0.004],
                    [0.002, 0.004, 0.025],
                ],
                dtype=float,
            ),
            numpy.array([-0.8, 0.9, 0.1], dtype=float),
        ),
    ],
)
def test_msd_latt_anisotropic(U, cartesian_vector):
    """Check Atom.msd_latt() anisotropic coordinate-invariance."""
    hexagonal = Lattice(1, 1, 1, 90, 90, 120)
    atom = Atom("C", [0, 0, 0], lattice=hexagonal, U=U)
    lattice_vector = hexagonal.fractional(cartesian_vector)
    actual = atom.msd_latt(lattice_vector)
    expected = pytest.approx(atom.msd_cart(cartesian_vector), rel=1e-13, abs=1e-13)
    assert actual == expected


@pytest.mark.parametrize(
    "uiso, cartesian_vector",
    [  # C1: isotropic displacement with msd is direction-independent in cartesian coordinates
        # Expected msd_cart equals Uisoequiv for any direction
        (0.0456, [0, 5, -2]),
    ],
)
def test_msd_cart_isotropic(uiso, cartesian_vector):
    """Check Atom.msd_cart()."""
    hexagonal = Lattice(1, 1, 1, 90, 90, 120)
    atom = Atom("C", [0, 0, 0], lattice=hexagonal, Uisoequiv=uiso)

    actual = atom.msd_cart(cartesian_vector)
    expected = pytest.approx(uiso, rel=0, abs=1e-15)
    assert actual == expected


@pytest.mark.parametrize(
    "U, cartesian_vector, scale",
    [  # C2: anisotropic displacement with msd_cart normalizes direction vector internally
        # Expected msd_cart(cartesian_vector) == msd_cart(scale * cartesian_vector)
        (
            numpy.array(
                [
                    [0.011, 0.001, 0.000],
                    [0.001, 0.019, 0.002],
                    [0.000, 0.002, 0.027],
                ],
                dtype=float,
            ),
            numpy.array([0.4, 1.1, -0.6], dtype=float),
            3.7,
        ),
        (
            numpy.array(
                [
                    [0.020, 0.003, -0.001],
                    [0.003, 0.014, 0.002],
                    [-0.001, 0.002, 0.009],
                ],
                dtype=float,
            ),
            numpy.array([2.0, -1.0, 0.5], dtype=float),
            0.25,
        ),
    ],
)
def test_msd_cart_anisotropic(U, cartesian_vector, scale):
    """Check Atom.msd_cart() anisotropic normalization invariance."""
    hexagonal = Lattice(1, 1, 1, 90, 90, 120)
    atom = Atom("C", [0, 0, 0], lattice=hexagonal, U=U)
    actual = atom.msd_cart(cartesian_vector)
    expected = pytest.approx(atom.msd_cart(scale * cartesian_vector), rel=1e-13, abs=1e-13)
    assert actual == expected


if __name__ == "__main__":
    unittest.main()

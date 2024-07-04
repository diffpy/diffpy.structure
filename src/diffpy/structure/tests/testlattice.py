#!/usr/bin/env python
##############################################################################
#
# diffpy.structure  by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2006 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################

"""Unit tests for Lattice class.
"""

import unittest

import numpy
import numpy.linalg as numalg

from diffpy.structure import Lattice, LatticeError

# ----------------------------------------------------------------------------


class TestLattice(unittest.TestCase):
    """test methods of Lattice class"""

    def setUp(self):
        self.lattice = Lattice()
        self.places = 12
        return

    def test___init__(self):
        """Check Lattice.__init__ processing of arguments."""
        self.assertRaises(ValueError, Lattice, self.lattice, c=4)
        self.assertRaises(ValueError, Lattice, base=self.lattice.base, baserot=self.lattice.baserot)
        self.assertRaises(ValueError, Lattice, 1, 2, 3)
        self.assertRaises(ValueError, Lattice, 1, 2, 3, 80, 90)
        L0 = self.lattice
        L0.setLatBase(L0.cartesian([[1, 1, 0], [0, 1, 1], [1, 0, 1]]))
        L1 = Lattice(L0)
        self.assertTrue(numpy.array_equal(L0.base, L1.base))
        L2 = Lattice(base=L0.base)
        self.assertTrue(numpy.array_equal(L0.base, L2.base))
        self.assertTrue(numpy.array_equal(L0.isotropicunit, L2.isotropicunit))
        L3 = Lattice(*L0.abcABG(), baserot=L0.baserot)
        self.assertTrue(numpy.allclose(L0.base, L3.base))
        self.assertTrue(numpy.allclose(L0.isotropicunit, L3.isotropicunit))
        return

    def test_setLatPar(self):
        """check calculation of standard unit cell vectors"""
        from math import cos, radians, sqrt

        from numpy import dot

        def norm(x):
            return sqrt(sum([xi**2 for xi in x]))

        def cosd(x):
            return cos(radians(x))

        self.lattice.setLatPar(1.0, 2.0, 3.0, 80, 100, 120)
        base = self.lattice.base
        self.assertAlmostEqual(1.0, norm(base[0]), self.places)
        self.assertAlmostEqual(2.0, norm(base[1]), self.places)
        self.assertAlmostEqual(3.0, norm(base[2]), self.places)
        self.assertAlmostEqual(cosd(80.0), dot(base[1], base[2]) / (2 * 3), self.places)
        self.assertAlmostEqual(cosd(100.0), dot(base[0], base[2]) / (1 * 3), self.places)
        self.assertAlmostEqual(cosd(120.0), dot(base[0], base[1]) / (1 * 2), self.places)
        return

    def test_latpar_properties(self):
        """check assignment to a, b, c, alpha, beta, gamma."""
        lat = self.lattice
        lat.a = 2
        lat.b = 4
        lat.c = 6
        lat.alpha = 80
        lat.beta = 100
        lat.gamma = 120
        lat1 = Lattice(2, 4, 6, 80, 100, 120)
        self.assertAlmostEqual(-0.5, lat.cg, self.places)
        self.assertTrue(numpy.array_equal(lat1.base, lat.base))
        return

    def test_readonly_properties(self):
        """Check that read-only properties are indeed such."""
        lat = self.lattice
        lat.b = 2
        lat.c = 6
        self.assertEqual(1.0, lat.unitvolume)
        self.assertRaises(AttributeError, setattr, lat, "unitvolume", 3.33)
        self.assertEqual(12, lat.volume)
        self.assertRaises(AttributeError, setattr, lat, "volume", 3.33)
        self.assertEqual(0.0, lat.ca)
        self.assertRaises(AttributeError, setattr, lat, "ca", 3.33)
        self.assertEqual(0.0, lat.cb)
        self.assertRaises(AttributeError, setattr, lat, "cb", 3.33)
        self.assertEqual(0.0, lat.cg)
        self.assertRaises(AttributeError, setattr, lat, "cg", 3.33)
        self.assertEqual(1.0, lat.sa)
        self.assertRaises(AttributeError, setattr, lat, "sa", 3.33)
        self.assertEqual(1.0, lat.sb)
        self.assertRaises(AttributeError, setattr, lat, "sb", 3.33)
        self.assertEqual(1.0, lat.sg)
        self.assertRaises(AttributeError, setattr, lat, "sg", 3.33)
        self.assertEqual(1.0, lat.ar)
        self.assertRaises(AttributeError, setattr, lat, "ar", 3.33)
        self.assertEqual(0.5, lat.br)
        self.assertRaises(AttributeError, setattr, lat, "br", 3.33)
        self.assertAlmostEqual(1.0 / 6, lat.cr, self.places)
        self.assertRaises(AttributeError, setattr, lat, "cr", 3.33)
        self.assertEqual(90.0, lat.alphar)
        self.assertRaises(AttributeError, setattr, lat, "alphar", 3.33)
        self.assertEqual(90.0, lat.betar)
        self.assertRaises(AttributeError, setattr, lat, "betar", 3.33)
        self.assertEqual(90.0, lat.gammar)
        self.assertRaises(AttributeError, setattr, lat, "gammar", 3.33)
        self.assertEqual(0.0, lat.car)
        self.assertRaises(AttributeError, setattr, lat, "car", 3.33)
        self.assertEqual(0.0, lat.cbr)
        self.assertRaises(AttributeError, setattr, lat, "cbr", 3.33)
        self.assertEqual(0.0, lat.cgr)
        self.assertRaises(AttributeError, setattr, lat, "cgr", 3.33)
        self.assertEqual(1.0, lat.sar)
        self.assertRaises(AttributeError, setattr, lat, "sar", 3.33)
        self.assertEqual(1.0, lat.sbr)
        self.assertRaises(AttributeError, setattr, lat, "sbr", 3.33)
        self.assertEqual(1.0, lat.sgr)
        self.assertRaises(AttributeError, setattr, lat, "sgr", 3.33)
        return

    def test_setLatBase(self):
        """check calculation of unit cell rotation"""
        base = numpy.array([[1.0, 1.0, 0.0], [0.0, 1.0, 1.0], [1.0, 0.0, 1.0]])
        self.lattice.setLatBase(base)
        self.assertAlmostEqual(self.lattice.a, numpy.sqrt(2.0), self.places)
        self.assertAlmostEqual(self.lattice.b, numpy.sqrt(2.0), self.places)
        self.assertAlmostEqual(self.lattice.c, numpy.sqrt(2.0), self.places)
        self.assertAlmostEqual(self.lattice.alpha, 60.0, self.places)
        self.assertAlmostEqual(self.lattice.beta, 60.0, self.places)
        self.assertAlmostEqual(self.lattice.gamma, 60.0, self.places)
        detR0 = numalg.det(self.lattice.baserot)
        self.assertAlmostEqual(detR0, 1.0, self.places)
        # try if rotation matrix works
        self.assertEqual(numpy.all(base == self.lattice.base), True)
        self.lattice.setLatPar(alpha=44, beta=66, gamma=88)
        self.assertNotEqual(numpy.all(base == self.lattice.base), True)
        self.lattice.setLatPar(alpha=60, beta=60, gamma=60)
        self.assertTrue(numpy.allclose(base[0], self.lattice.base[0]))
        self.assertTrue(numpy.allclose(base[1], self.lattice.base[1]))
        self.assertTrue(numpy.allclose(base[2], self.lattice.base[2]))
        # try base checking
        self.assertRaises(LatticeError, self.lattice.setLatBase, [[1, 0, 0], [1, 0, 0], [0, 0, 1]])
        self.assertRaises(LatticeError, self.lattice.setLatBase, [[1, 0, 0], [0, 0, 1], [0, 1, 0]])
        return

    def test_reciprocal(self):
        """check calculation of reciprocal lattice."""
        r1 = self.lattice.reciprocal()
        self.assertEqual((1, 1, 1, 90, 90, 90), r1.abcABG())
        L2 = Lattice(2, 4, 8, 90, 90, 90)
        r2 = L2.reciprocal()
        self.assertEqual((0.5, 0.25, 0.125, 90, 90, 90), r2.abcABG())
        rr2 = r2.reciprocal()
        self.assertTrue(numpy.array_equal(L2.base, rr2.base))
        return

    def test_dot(self):
        """check dot product of lattice vectors."""
        L = self.lattice
        L.setLatPar(gamma=120)
        self.assertAlmostEqual(-0.5, L.dot([1, 0, 0], [0, 1, 0]), self.places)
        va5 = numpy.tile([1.0, 0.0, 0.0], (5, 1))
        vb5 = numpy.tile([0.0, 1.0, 0.0], (5, 1))
        self.assertTrue(numpy.array_equal(5 * [-0.5], L.dot(va5, vb5)))
        self.assertTrue(numpy.array_equal(5 * [-0.5], L.dot(va5[0], vb5)))
        self.assertTrue(numpy.array_equal(5 * [-0.5], L.dot(va5, vb5[0])))
        return

    def test_norm(self):
        """check norm of a lattice vector."""
        self.assertEqual(1, self.lattice.norm([1, 0, 0]))
        u = numpy.array([[3, 4, 0], [1, 1, 1]])
        self.assertTrue(numpy.allclose([5, 3**0.5], self.lattice.norm(u)))
        self.lattice.setLatPar(gamma=120)
        self.assertAlmostEqual(1, self.lattice.norm([1, 1, 0]), self.places)
        return

    def test_rnorm(self):
        """check norm of a reciprocal vector."""
        L = self.lattice
        L.setLatPar(1, 1.5, 2.3, 80, 95, 115)
        r = L.reciprocal()
        hkl = [0.5, 0.3, 0.2]
        self.assertAlmostEqual(r.norm(hkl), L.rnorm(hkl), self.places)
        hkl5 = numpy.tile(hkl, (5, 1))
        self.assertTrue(numpy.allclose(5 * [r.norm(hkl)], L.rnorm(hkl5)))
        return

    def test_dist(self):
        """check dist function for distance between lattice points."""
        L = self.lattice
        L.setLatPar(1, 1.5, 2.3, 80, 95, 115)
        u = [0.1, 0.3, 0.7]
        v = [0.3, 0.7, 0.7]
        d0 = numalg.norm(L.cartesian(numpy.array(u) - v))
        self.assertAlmostEqual(d0, L.dist(u, v), self.places)
        self.assertAlmostEqual(d0, L.dist(v, u), self.places)
        u5 = numpy.tile(u, (5, 1))
        v5 = numpy.tile(v, (5, 1))
        self.assertTrue(numpy.allclose(5 * [d0], L.dist(u, v5)))
        self.assertTrue(numpy.allclose(5 * [d0], L.dist(u5, v)))
        self.assertTrue(numpy.allclose(5 * [d0], L.dist(v5, u5)))
        return

    def test_angle(self):
        """check angle calculation between lattice vectors."""
        from math import acos, degrees

        L = self.lattice
        L.setLatPar(1, 1.5, 2.3, 80, 95, 115)
        u = [0.1, 0.3, 0.7]
        v = [0.3, 0.7, 0.7]
        uc = L.cartesian(u)
        vc = L.cartesian(v)
        a0 = degrees(acos(numpy.dot(uc, vc) / (numalg.norm(uc) * numalg.norm(vc))))
        self.assertAlmostEqual(a0, L.angle(u, v), self.places)
        self.assertAlmostEqual(a0, L.angle(v, u), self.places)
        u5 = numpy.tile(u, (5, 1))
        v5 = numpy.tile(v, (5, 1))
        self.assertTrue(numpy.allclose(5 * [a0], L.angle(u, v5)))
        self.assertTrue(numpy.allclose(5 * [a0], L.angle(u5, v)))
        self.assertTrue(numpy.allclose(5 * [a0], L.angle(v5, u5)))
        return

    def test_repr(self):
        """check string representation of this lattice"""
        r = repr(self.lattice)
        self.assertEqual(r, "Lattice()")
        self.lattice.setLatPar(1, 2, 3, 10, 20, 30)
        r = repr(self.lattice)
        r0 = "Lattice(a=1, b=2, c=3, alpha=10, beta=20, gamma=30)"
        self.assertEqual(r, r0)
        base = [[1.0, 1.0, 0.0], [0.0, 2.0, 2.0], [3.0, 0.0, 3.0]]
        self.lattice.setLatBase(base)
        r = repr(self.lattice)
        self.assertEqual(r, "Lattice(base=%r)" % self.lattice.base)


# End of class TestLattice

# ----------------------------------------------------------------------------

if __name__ == "__main__":
    unittest.main()

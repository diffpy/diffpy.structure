##############################################################################
#
# Structure         by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2006 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################

"""Unit tests for Lattice class.
"""

__id__ = "$Id$"

import unittest

from diffpy.Structure import Lattice, LatticeError

##############################################################################
class TestLattice(unittest.TestCase):
    """test methods of Lattice class"""

    def setUp(self):
        self.lattice = Lattice()
        self.places = 12
        return

    def assertListAlmostEqual(self, l1, l2, places=None):
        """wrapper for list comparison"""
        if places is None: places = self.places
        self.assertEqual(len(l1), len(l2))
        for i in range(len(l1)):
            self.assertAlmostEqual(l1[i], l2[i], places)

    def test_setLatPar(self):
        """check calculation of standard unit cell vectors"""
        from numpy import dot
        from math import radians, sqrt, cos, sin
        norm = lambda x : sqrt(sum([xi**2 for xi in x]))
        cosd = lambda x : cos(radians(x))
        sind = lambda x : sin(radians(x))
        self.lattice.setLatPar(1.0, 2.0, 3.0, 80, 100, 120)
        base = self.lattice.base
        self.assertAlmostEqual(1.0, norm(base[0]), self.places)
        self.assertAlmostEqual(2.0, norm(base[1]), self.places)
        self.assertAlmostEqual(3.0, norm(base[2]), self.places)
        self.assertAlmostEqual(cosd(80.0),
                dot(base[1],base[2])/(2*3), self.places)
        self.assertAlmostEqual(cosd(100.0),
                dot(base[0],base[2])/(1*3), self.places)
        self.assertAlmostEqual(cosd(120.0),
                dot(base[0],base[1])/(1*2), self.places)

    def test_setLatBase(self):
        """check calculation of unit cell rotation"""
        import numpy
        import numpy.linalg as numalg
        base = numpy.array([[ 1.0,  1.0,  0.0],
                          [ 0.0,  1.0,  1.0],
                          [ 1.0,  0.0,  1.0]])
        self.lattice.setLatBase(base)
        self.assertAlmostEqual(self.lattice.a, numpy.sqrt(2.0), self.places)
        self.assertAlmostEqual(self.lattice.b, numpy.sqrt(2.0), self.places)
        self.assertAlmostEqual(self.lattice.c, numpy.sqrt(2.0), self.places)
        self.assertAlmostEqual(self.lattice.alpha, 60.0, self.places)
        self.assertAlmostEqual(self.lattice.beta,  60.0, self.places)
        self.assertAlmostEqual(self.lattice.gamma, 60.0, self.places)
        detR0 = numalg.det(self.lattice.baserot)
        self.assertAlmostEqual(detR0, 1.0, self.places)
        # try if rotation matrix works
        self.assertEqual(numpy.all(base == self.lattice.base), True)
        self.lattice.setLatPar(alpha=44, beta=66, gamma=88)
        self.assertNotEqual(numpy.all(base == self.lattice.base), True)
        self.lattice.setLatPar(alpha=60, beta=60, gamma=60)
        self.assertListAlmostEqual(base[0], self.lattice.base[0])
        self.assertListAlmostEqual(base[1], self.lattice.base[1])
        self.assertListAlmostEqual(base[2], self.lattice.base[2])
        # try base checking
        self.assertRaises(LatticeError, self.lattice.setLatBase,
                [[1, 0, 0], [1,0,0], [0,0,1]])
        self.assertRaises(LatticeError, self.lattice.setLatBase,
                [[1, 0, 0], [0,0,1], [0,1,0]])
        return

    def test_repr(self):
        """check string representation of this lattice"""
        r = repr(self.lattice)
        self.assertEqual(r, "Lattice()")
        self.lattice.setLatPar(1, 2, 3, 10, 20, 30)
        r = repr(self.lattice)
        r0 = "Lattice(a=1, b=2, c=3, alpha=10, beta=20, gamma=30)"
        self.assertEqual(r, r0)
        base = [[ 1.0,  1.0,  0.0],
                [ 0.0,  2.0,  2.0],
                [ 3.0,  0.0,  3.0]]
        self.lattice.setLatBase(base)
        r = repr(self.lattice)
        self.assertEqual(r, "Lattice(base=%r)" % self.lattice.base)

# End of TestLattice

if __name__ == '__main__':
    unittest.main()

# End of file

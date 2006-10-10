#!/usr/bin/env python

"""Unit tests for SymmetryUtilities.py
"""

# version
__id__ = '$Id$'

import os
import sys
import unittest

import Structure.SpaceGroups
from Structure.SymmetryUtilities import *

##############################################################################
class TestRoutines(unittest.TestCase):

    def setUp(self):
        return

    def tearDown(self):
        return

    def test_isSpaceGroupLatPar(self):
        """check isSpaceGroupLatPar()
        """
        triclinic = Structure.SpaceGroups.GetSpaceGroup("P1")
        monoclinic = Structure.SpaceGroups.GetSpaceGroup("P2")
        orthorhombic = Structure.SpaceGroups.GetSpaceGroup("P222")
        tetragonal = Structure.SpaceGroups.GetSpaceGroup("P4")
        trigonal = Structure.SpaceGroups.GetSpaceGroup("P3") 
        hexagonal = Structure.SpaceGroups.GetSpaceGroup("P6") 
        cubic = Structure.SpaceGroups.GetSpaceGroup("P23")
        self.failUnless(isSpaceGroupLatPar(triclinic, 1, 2, 3, 40, 50, 60))
        self.failIf(isSpaceGroupLatPar(monoclinic, 1, 2, 3, 40, 50, 60))
        self.failUnless(isSpaceGroupLatPar(monoclinic, 1, 2, 3, 90, 50, 90))
        self.failIf(isSpaceGroupLatPar(orthorhombic, 1, 2, 3, 90, 50, 90))
        self.failUnless(isSpaceGroupLatPar(orthorhombic, 1, 2, 3, 90, 90, 90))
        self.failIf(isSpaceGroupLatPar(tetragonal, 1, 2, 3, 90, 90, 90))
        self.failUnless(isSpaceGroupLatPar(tetragonal, 2, 2, 3, 90, 90, 90))
        self.failIf(isSpaceGroupLatPar(trigonal, 2, 2, 3, 90, 90, 90))
        self.failUnless(isSpaceGroupLatPar(trigonal, 2, 2, 2, 80, 80, 80))
        self.failIf(isSpaceGroupLatPar(hexagonal, 2, 2, 2, 80, 80, 80))
        self.failUnless(isSpaceGroupLatPar(hexagonal, 2, 2, 3, 90, 90, 120))
        self.failIf(isSpaceGroupLatPar(cubic, 2, 2, 3, 90, 90, 120))
        self.failUnless(isSpaceGroupLatPar(cubic, 3, 3, 3, 90, 90, 90))
        return

    def test_expandPosition(self):
        """check expandPosition()
        """
        # ok again Ni example
        fcc = Structure.SpaceGroups.GetSpaceGroup(225)
        pos,pops,pmult = expandPosition(fcc, [0,0,0])
        self.failUnless(numpy.all(pos[0] == 0.0))
        self.assertEqual(4, len(pos))
        self.assertEqual(192, sum([len(l) for l in pops]))
        self.assertEqual(4, pmult)
        return

# End of class TestRoutines

##############################################################################
class TestPosition2Tuple(unittest.TestCase):

    def setUp(self):
        self.eps = 1.0e-4
        self.pos2tuple = Position2Tuple(self.eps)
        return

    def tearDown(self):
        del self.pos2tuple
        return

    def test___init__(self):
        """check Position2Tuple.__init__()
        """
        self.assertNotEqual(0.0, self.pos2tuple.eps)
        self.pos2tuple = Position2Tuple(1.0/sys.maxint/2)
        self.assertEqual(0.0, self.pos2tuple.eps)
        return

    def test___call__(self):
        """check Position2Tuple.__call__()
        """
        pos2tuple = self.pos2tuple
        import numpy
        positions = numpy.zeros((100,3), dtype=numpy.float64)
        positions[:,0] = numpy.arange(100)/100.0*pos2tuple.eps + 0.1
        positions = positions - numpy.floor(positions)
        # pos2tuple should generate at most 2 distinct tuples
        alltuples = dict.fromkeys([pos2tuple(xyz) for xyz in positions])
        self.failIf(len(alltuples) > 2)
        return

# End of class TestPosition2Tuple

if __name__ == '__main__':
    unittest.main()

# End of file

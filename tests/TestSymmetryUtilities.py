########################################################################
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
########################################################################

"""Unit tests for SymmetryUtilities.py
"""

# version
__id__ = '$Id$'

import os
import sys
import unittest

from diffpy.Structure.SpaceGroups import GetSpaceGroup
from diffpy.Structure.SymmetryUtilities import *

##############################################################################
class TestRoutines(unittest.TestCase):

    def setUp(self):
        return

    def tearDown(self):
        return

    def test_isSpaceGroupLatPar(self):
        """check isSpaceGroupLatPar()
        """
        triclinic = GetSpaceGroup("P1")
        monoclinic = GetSpaceGroup("P2")
        orthorhombic = GetSpaceGroup("P222")
        tetragonal = GetSpaceGroup("P4")
        trigonal = GetSpaceGroup("P3") 
        hexagonal = GetSpaceGroup("P6") 
        cubic = GetSpaceGroup("P23")
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
        fcc = GetSpaceGroup(225)
        pos,pops,pmult = expandPosition(fcc, [0,0,0])
        self.failUnless(numpy.all(pos[0] == 0.0))
        self.assertEqual(4, len(pos))
        self.assertEqual(192, sum([len(l) for l in pops]))
        self.assertEqual(4, pmult)
        return

    def test_pruneFormulaDictionary(self):
        """check pruneFormulaDictionary()
        """
        fmdict = {"x" : "3*y-0.17", "y" : '0', "z" : "0.13"}
        pruned = pruneFormulaDictionary(fmdict)
        self.assertEqual({"x" : "3*y-0.17"}, pruned)
        return

    def test_isconstantFormula(self):
        """check isconstantFormula()
        """
        self.failIf(isconstantFormula('x-y+z'))
        self.failUnless(isconstantFormula('6.023e23'))
        self.failUnless(isconstantFormula('22/7'))
        self.failUnless(isconstantFormula('- 22/7'))
        self.failUnless(isconstantFormula('+13/ 9'))
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
        positions = numpy.zeros((100,3), dtype=float)
        positions[:,0] = numpy.arange(100)/100.0*pos2tuple.eps + 0.1
        positions = positions - numpy.floor(positions)
        # pos2tuple should generate at most 2 distinct tuples
        alltuples = dict.fromkeys([pos2tuple(xyz) for xyz in positions])
        self.failIf(len(alltuples) > 2)
        return

# End of class TestPosition2Tuple

##############################################################################
class TestGeneratorSite(unittest.TestCase):

    generators = {}

    def setUp(self):
        x, y, z = 0.07, 0.11, 0.13
        self.x, self.y, self.z = x, y, z 
        if TestGeneratorSite.generators:
            self.__dict__.update(TestGeneratorSite.generators)
            return
        sg117 = GetSpaceGroup(117)
        sg143 = GetSpaceGroup(143)
        sg227 = GetSpaceGroup(227)
        g117c = GeneratorSite(sg117, [0, 0.5, 0])
        g117h = GeneratorSite(sg117, [x, x+0.5, 0.5])
        g143a = GeneratorSite(sg143, [0, 0, z])
        g143b = GeneratorSite(sg143, [1./3, 2./3, z])
        g143c = GeneratorSite(sg143, [2./3, 1./3, z])
        g143d = GeneratorSite(sg143, [x, y, z])
        g227a = GeneratorSite(sg227, [0, 0, 0])
        g227c = GeneratorSite(sg227, 3*[1./8])
        g227oa = GeneratorSite(sg227, 3*[1./8], sgoffset=3*[1./8])
        g227oc = GeneratorSite(sg227, [0, 0, 0], sgoffset=3*[1./8])
        TestGeneratorSite.generators = {
                'g117c' : g117c,  'g117h' : g117h, 
                'g143a' : g143a,  'g143b' : g143b,
                'g143c' : g143c,  'g143d' : g143d,
                'g227a' : g227a,  'g227c' : g227c,
                'g227oa' : g227oa,  'g227oc' : g227oc
        }
        self.__dict__.update(TestGeneratorSite.generators)
        return

    def tearDown(self):
        return

    def test___init__(self):
        """check GeneratorSite.__init__()
        """
        # check multiplicities
        self.assertEqual(2,  self.g117c.multiplicity)
        self.assertEqual(4,  self.g117h.multiplicity)
        self.assertEqual(1,  self.g143a.multiplicity)
        self.assertEqual(1,  self.g143b.multiplicity)
        self.assertEqual(1,  self.g143c.multiplicity)
        self.assertEqual(3,  self.g143d.multiplicity)
        self.assertEqual(8,  self.g227a.multiplicity)
        self.assertEqual(16, self.g227c.multiplicity)
        self.assertEqual(8,  self.g227oa.multiplicity)
        self.assertEqual(16, self.g227oc.multiplicity)
        return

    def test_positionFormula(self):
        """check GeneratorSite.positionFormula()
        """
        import re
        # 117c
        self.assertEqual([], self.g117c.pparameters)
        self.assertEqual([("x", self.x)], self.g117h.pparameters)
        # 143c
        pfm143c = self.g143c.positionFormula(self.g143c.xyz)
        places = 12
        self.assertEqual("+2/3", pfm143c["x"])
        self.assertEqual("+1/3", pfm143c["y"])
        self.assertEqual("z", pfm143c["z"])
        # 143d
        x, y, z = self.x, self.y, self.z
        pfm143d = self.g143d.positionFormula([-x+y, -x, z])
        self.assertEqual("-x+y", pfm143d["x"].replace(' ',''))
        self.assertEqual("-x", pfm143d["y"].replace(' ',''))
        self.failUnless(re.match("[+]?z", pfm143d["z"].strip()))
        # 227a
        self.assertEqual([], self.g227a.pparameters)
        self.assertEqual([], self.g227oa.pparameters)
        # 227c
        self.assertEqual([], self.g227c.pparameters)
        self.assertEqual([], self.g227oc.pparameters)
        return

    def test_UFormula(self):
        """check GeneratorSite.UFormula()
        """
        # Ref: Willis and Pryor, Thermal Vibrations in Crystallography,
        # Cambridge University Press 1975, p. 104-110
        smbl = ('A', 'B', 'C', 'D', 'E', 'F')
        norule = { 'U11':'A', 'U22':'B', 'U33':'C',
                   'U12':'D', 'U13':'E', 'U23':'F' }
        rule05 = { 'U11':'A', 'U22':'A', 'U33':'C',
                   'U12':'D', 'U13':'0', 'U23':'0' }
        rule07 = { 'U11':'A', 'U22':'A', 'U33':'C',
                   'U12':'D', 'U13':'E', 'U23':'-E' }
        rule16 = { 'U11':'A', 'U22':'A', 'U33':'C',
                   'U12':'A', 'U13':'0', 'U23':'0' }
        rule17 = { 'U11':'A', 'U22':'A', 'U33':'A',
                   'U12':'0', 'U13':'0', 'U23':'0' }
        rule18 = { 'U11':'A', 'U22':'A', 'U33':'A',
                   'U12':'D', 'U13':'D', 'U23':'D' }
        ufm = self.g117c.UFormula(self.g117c.xyz, smbl)
        self.assertEqual(rule05, ufm)
        ufm = self.g117h.UFormula(self.g117h.xyz, smbl)
        self.assertEqual(rule07, ufm)
        ufm = self.g143a.UFormula(self.g143a.xyz, smbl)
        self.assertEqual(rule16, ufm)
        ufm = self.g143b.UFormula(self.g143b.xyz, smbl)
        self.assertEqual(rule16, ufm)
        ufm = self.g143c.UFormula(self.g143c.xyz, smbl)
        self.assertEqual(rule16, ufm)
        ufm = self.g143d.UFormula(self.g143d.xyz, smbl)
        self.assertEqual(norule, ufm)
        ufm = self.g227a.UFormula(self.g227a.xyz, smbl)
        self.assertEqual(rule17, ufm)
        ufm = self.g227c.UFormula(self.g227c.xyz, smbl)
        self.assertEqual(rule18, ufm)
        ufm = self.g227oa.UFormula(self.g227oa.xyz, smbl)
        self.assertEqual(rule17, ufm)
        ufm = self.g227oc.UFormula(self.g227oc.xyz, smbl)
        self.assertEqual(rule18, ufm)
        return

    def test_eqIndex(self):
        """check GeneratorSite.eqIndex()
        """
        self.assertEqual(13, self.g227oc.eqIndex(self.g227oc.eqxyz[13]))
        return

# End of class TestGeneratorSite

if __name__ == '__main__':
    unittest.main()

# End of file

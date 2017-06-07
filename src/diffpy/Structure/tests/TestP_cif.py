#!/usr/bin/env python
##############################################################################
#
# diffpy.Structure  by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2007 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################

"""Unit tests for diffpy.Structure.Parsers.P_cif module
"""

import unittest
import numpy

from diffpy.Structure.tests.testutils import datafile
from diffpy.Structure.Parsers.P_cif import P_cif, leading_float, getSymOp
from diffpy.Structure.Parsers import getParser
from diffpy.Structure import Structure
from diffpy.Structure import StructureFormatError

##############################################################################
class TestRoutines(unittest.TestCase):

    def setUp(self):
        return

    def tearDown(self):
        return

    def test_leading_float(self):
        """check leading_float()
        """
        self.assertEqual(0.37, leading_float('0.37(3)'))
        self.assertEqual(0.37, leading_float('0.37ab\ncd'))
        self.assertRaises(ValueError, leading_float, 'q1')
        return

    def test_getSymOp(self):
        """check getSymOp()
        """
        from diffpy.Structure.SpaceGroups import SymOp
        from diffpy.Structure.SpaceGroups import Rot_X_mY_Z, Tr_0_12_12
        op = getSymOp('x,1/2-y,1/2+z')
        op_std = SymOp(Rot_X_mY_Z, Tr_0_12_12)
        self.assertEqual(str(op_std), str(op))
        from diffpy.Structure.SpaceGroups import Rot_mX_mXY_Z, Tr_0_0_12
        op1 = getSymOp('-x,-x+y,1/2+z')
        op1_std = SymOp(Rot_mX_mXY_Z, Tr_0_0_12)
        self.assertEqual(str(op1_std), str(op1))
        return

# End of class TestRoutines

##############################################################################
class TestP_cif(unittest.TestCase):

    goodciffile = datafile('PbTe.cif')
    badciffile = datafile('LiCl-bad.cif')
    graphiteciffile = datafile('graphite.cif')
    cdsebulkpdffitfile = datafile('CdSe_bulk.stru')
    teiciffile = datafile('TeI.cif')
    places = 6

    def setUp(self):
        self.ptest = P_cif()
        self.pfile = P_cif()
        return

    def tearDown(self):
        return

    def test_parse(self):
        """check P_cif.parse()
        """
        with open(self.goodciffile) as fp1:
            sgood = fp1.read()
        with open(self.badciffile) as fp2:
            sbad = fp2.read()
        pfile, ptest = self.pfile, self.ptest
        stru_check = pfile.parseFile(self.goodciffile)
        stru = ptest.parse(sgood)
        self.assertEqual(str(stru_check), str(stru))
        self.assertEqual(str(stru_check.lattice), str(stru.lattice))
        self.assertEqual(pfile.spacegroup.short_name,
            ptest.spacegroup.short_name)
        ptestb = P_cif()
        self.assertRaises(StructureFormatError,
            ptestb.parse, sbad)
        return

    def test_parseLines(self):
        """check P_cif.parseLines()
        """
        with open(self.goodciffile) as fp1:
            goodlines = fp1.readlines()
        with open(self.badciffile) as fp2:
            badlines = fp2.readlines()
        pfile, ptest = self.pfile, self.ptest
        stru_check = pfile.parseFile(self.goodciffile)
        stru = ptest.parseLines(goodlines)
        self.assertEqual(str(stru_check), str(stru))
        self.assertEqual(str(stru_check.lattice), str(stru.lattice))
        self.assertEqual(pfile.spacegroup.short_name,
            ptest.spacegroup.short_name)
        ptest2 = P_cif()
        self.assertRaises(StructureFormatError,
                ptest2.parseLines, badlines)
        return

    def test_parseFile(self):
        """check P_cif.parseFile()
        """
        # goodciffile
        stru = self.pfile.parseFile(self.goodciffile)
        self.assertEqual(8, len(stru))
        self.assertEqual(6.461, stru.lattice.a)
        self.assertEqual(6.461, stru.lattice.b)
        self.assertEqual(6.461, stru.lattice.c)
        self.assertEqual(90.0, stru.lattice.alpha)
        self.assertEqual(90.0, stru.lattice.beta)
        self.assertEqual(90.0, stru.lattice.gamma)
        self.assertEqual('Fm-3m', self.pfile.spacegroup.short_name)
        a0 = stru[0]
        self.assertEqual(0.5, a0.x)
        self.assertEqual(0.5, a0.y)
        self.assertEqual(0.5, a0.z)
        self.assertEqual(False, a0.anisotropy)
        self.assertEqual(1.0, a0.occupancy)
        self.assertEqual(0.0225566, a0.Uisoequiv)
        # badciffile
        pfile2 = P_cif()
        self.assertRaises(StructureFormatError,
                pfile2.parseFile, self.badciffile)
        # graphite
        pgraphite = P_cif()
        graphite = pgraphite.parseFile(self.graphiteciffile)
        self.assertEqual(4, len(graphite))
        c1 = graphite[0]
        self.assertEqual(str, type(c1.element))
        self.assertEqual('C', c1.element)
        self.assertEqual(str, type(c1.label))
        self.assertEqual('C1', c1.label)
        # filename with unicode encoding
        ugraphite = P_cif().parseFile(unicode(self.graphiteciffile))
        self.assertEqual(4, len(ugraphite))
        # File with full space group name
        ptei = P_cif().parseFile(self.teiciffile)
        self.assertEqual(16 , len(ptei))
        return

#   def test__parseCifBlock(self):
#       """check P_cif._parseCifBlock()
#       """
#       return
#
#   def test__parse_lattice(self):
#       """check P_cif._parse_lattice()
#       """
#       return
#
#   def test__parse_atom_site_label(self):
#       """check P_cif._parse_atom_site_label()
#       """
#       return
#
#   def test__parse_atom_site_aniso_label(self):
#       """check P_cif._parse_atom_site_aniso_label()
#       """
#       return
#
#   def test__parse_space_group_symop_operation_xyz(self):
#       """check P_cif._parse_space_group_symop_operation_xyz()
#       """
#       return
#
#   def test__expandAsymmetricUnit(self):
#       """check P_cif._expandAsymmetricUnit()
#       """
#       return
#
#   def test_toLines(self):
#       """check P_cif.toLines()
#       """
#       return
#
#   def test_tostring(self):
#       """check P_cif.tostring()
#       """
#       return

    def test_write_and_read(self):
        """high-level check of P_cif.tostring()
        """
        # high-level check
        stru_check = Structure()
        stru_check.read(self.cdsebulkpdffitfile)
        s_s = stru_check.writeStr('cif')
        stru = Structure()
        stru.readStr(s_s, 'cif')
        self.assertAlmostEqual(4.2352, stru.lattice.a, self.places)
        self.assertAlmostEqual(4.2352, stru.lattice.b, self.places)
        self.assertAlmostEqual(6.90603, stru.lattice.c, self.places)
        self.assertEqual(4, len(stru))
        a0 = stru[0]
        self.assertEqual('Cd', a0.element)
        self.assertTrue(numpy.allclose([0.3334, 0.6667, 0.0], a0.xyz))
        self.assertTrue(a0.anisotropy)
        self.assertAlmostEqual(0.01303, a0.U[0,0])
        self.assertAlmostEqual(0.01303, a0.U[1,1])
        self.assertAlmostEqual(0.01402, a0.U[2,2])
        a3 = stru[3]
        self.assertEqual('Se', a3.element)
        self.assertTrue(numpy.allclose([0.6666, 0.333300, 0.87667], a3.xyz))
        self.assertAlmostEqual(0.015673, a3.U[0,0])
        self.assertAlmostEqual(0.015673, a3.U[1,1])
        self.assertAlmostEqual(0.046164, a3.U[2,2])
        return

    def test_eps(self):
        """Test the P_cif.eps coordinates resolution.
        """
        pcif = P_cif()
        pcif.eps = 1e-8
        grph = pcif.parseFile(self.graphiteciffile)
        self.assertEqual(8, len(grph))
        self.assertTrue(all(a.label.startswith('C1') for a in grph[:2]))
        self.assertTrue(all(a.label.startswith('C2') for a in grph[2:]))
        pcif2 = P_cif()
        pcif2.eps = 1e-3
        grph2 = pcif2.parseFile(self.graphiteciffile)
        self.assertEqual(4, len(grph2))
        return

    def test_getParser(self):
        """Test passing of eps keyword argument by getParser function.
        """
        pcif = getParser('cif', eps=1e-6)
        grph = pcif.parseFile(self.graphiteciffile)
        self.assertEqual(8, len(grph))
        self.assertTrue(all(a.label.startswith('C1') for a in grph[:2]))
        self.assertTrue(all(a.label.startswith('C2') for a in grph[2:]))
        pcif2 = getParser('cif')
        grph2 = pcif2.parseFile(self.graphiteciffile)
        self.assertEqual(4, len(grph2))
        return

# End of class TestP_cif

if __name__ == '__main__':
    unittest.main()

# End of file

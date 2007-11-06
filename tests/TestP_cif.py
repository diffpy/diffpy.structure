#!/usr/bin/env python

"""Unit tests for diffpy.Structure.Parsers.P_cif module
"""

# version
__id__ = '$Id$'

import os
import unittest

# useful variables
thisfile = locals().get('__file__', 'TestP_cif.py')
tests_dir = os.path.dirname(os.path.abspath(thisfile))
testdata_dir = os.path.join(tests_dir, 'testdata')

from diffpy.Structure.Parsers.P_cif import *
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
        from diffpy.Structure.SpaceGroups import SymOp, Rot_X_mY_Z, Tr_0_12_12
        op = getSymOp('x,1/2-y,1/2+z')
        op_std = SymOp(Rot_X_mY_Z, Tr_0_12_12)
        self.assertEqual(str(op_std), str(op))
        return

# End of class TestRoutines

##############################################################################
class TestP_cif(unittest.TestCase):

    def setUp(self):
        self.ptest = P_cif()
        self.pfile = P_cif()
        self.goodciffile = os.path.join(testdata_dir, 'PbTe.cif')
        self.badciffile = os.path.join(testdata_dir, 'LiCl-bad.cif')
        return

    def tearDown(self):
        return

    def test_parse(self):
        """check P_cif.parse()
        """
        sgood = open(self.goodciffile).read()
        sbad = open(self.badciffile).read()
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
        goodlines = open(self.goodciffile).readlines()
        badlines = open(self.badciffile).readlines()
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
        self.assertEqual(0.5, a0.xyz[0])
        self.assertEqual(0.5, a0.xyz[1])
        self.assertEqual(0.5, a0.xyz[2])
        self.assertEqual(False, a0.anisotropy)
        self.assertEqual(1.0, a0.occ)
        self.assertEqual(0.0225566, a0.Uisoequiv)
        pfile2 = P_cif()
        self.assertRaises(StructureFormatError,
                pfile2.parseFile, self.badciffile)
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

# End of class TestP_cif

if __name__ == '__main__':
    unittest.main()

# End of file

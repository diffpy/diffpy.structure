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

"""Unit tests for diffpy.structure.parsers.p_pdffit module
"""

import re
import unittest

import numpy
import pytest

from diffpy.structure import Structure
from diffpy.structure.structureerrors import StructureFormatError

# ----------------------------------------------------------------------------


class TestP_pdffit(unittest.TestCase):
    """test Parser for PDFFit file format"""

    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    def setUp(self):
        self.stru = Structure()
        self.format = "pdffit"
        self.places = 8

    def test_read_pdffit_ZnSb(self):
        """check reading of ZnSb pdffit structure file"""
        stru = self.stru
        stru.read(self.datafile("ZnSb_RT_Q28X_VM_20_fxiso.rstr"), self.format)
        f_title = "Cell structure file of Zn4Sb3 #167 interstitial"
        self.assertEqual(stru.title, f_title)
        self.assertAlmostEqual(stru.pdffit["scale"], 0.826145)
        self.assertAlmostEqual(stru.pdffit["delta2"], 4.687951)
        self.assertAlmostEqual(stru.pdffit["delta1"], 0.01)
        self.assertAlmostEqual(stru.pdffit["sratio"], 1.02)
        self.assertAlmostEqual(stru.pdffit["rcut"], 0.03)
        self.assertEqual(stru.pdffit["spcgr"], "R-3c")
        s_lat = [
            stru.lattice.a,
            stru.lattice.b,
            stru.lattice.c,
            stru.lattice.alpha,
            stru.lattice.beta,
            stru.lattice.gamma,
        ]
        f_lat = [12.309436, 12.309436, 12.392839, 90.0, 90.0, 120.0]
        self.assertTrue(numpy.allclose(s_lat, f_lat))
        s_dcell = stru.pdffit["dcell"]
        f_dcell = [0.000008, 0.000008, 0.000013, 0.0, 0.0, 0.0]
        self.assertTrue(numpy.allclose(s_dcell, f_dcell))
        self.assertEqual(stru.pdffit["ncell"], [1, 1, 1, 66])
        s_els = [a.element for a in stru]
        self.assertEqual(s_els, 36 * ["Zn"] + 30 * ["Sb"])
        a0 = stru[0]
        s_xyz = a0.xyz
        f_xyz = [0.09094387, 0.24639539, 0.40080261]
        s_o = a0.occupancy
        f_o = 0.9
        s_sigxyz = a0.sigxyz
        f_sigxyz = [0.00000079, 0.00000076, 0.00000064]
        s_sigo = a0.sigo
        f_sigo = 0.0
        s_U = [a0.U[i][i] for i in range(3)]
        f_U = 3 * [0.01]
        self.assertTrue(numpy.allclose(s_xyz, f_xyz))
        self.assertTrue(numpy.allclose(s_sigxyz, f_sigxyz))
        self.assertTrue(numpy.allclose(s_U, f_U))
        self.assertAlmostEqual(s_o, f_o)
        self.assertAlmostEqual(s_sigo, f_sigo)

    def test_read_pdffit_Ni(self):
        """check reading of Ni pdffit structure file"""
        stru = self.stru
        stru.read(self.datafile("Ni.stru"), self.format)
        f_title = "structure Ni  FCC"
        self.assertEqual(stru.title, f_title)
        self.assertEqual(stru.pdffit["spcgr"], "Fm-3m")
        s_lat = [
            stru.lattice.a,
            stru.lattice.b,
            stru.lattice.c,
            stru.lattice.alpha,
            stru.lattice.beta,
            stru.lattice.gamma,
        ]
        f_lat = [3.52, 3.52, 3.52, 90.0, 90.0, 90.0]
        for i in range(len(s_lat)):
            self.assertAlmostEqual(s_lat[i], f_lat[i])
        self.assertEqual(stru.pdffit["ncell"], [1, 1, 1, 4])
        s_els = [a.element for a in stru]
        self.assertEqual(s_els, 4 * ["Ni"])
        a0 = stru[0]
        s_xyz = a0.xyz
        f_xyz = [0.0, 0.0, 0.0]
        s_o = a0.occupancy
        f_o = 1.0
        s_U = [a0.U[i][i] for i in range(3)]
        f_U = 3 * [0.00126651]
        for i in range(3):
            self.assertAlmostEqual(s_xyz[i], f_xyz[i])
            self.assertAlmostEqual(s_U[i], f_U[i])
        self.assertAlmostEqual(s_o, f_o)

    def test_read_pdffit_Ni_prim123(self):
        """check reading of Ni_prim supercell 1x2x3"""
        stru = self.stru
        stru.read(self.datafile("Ni_prim123.stru"), self.format)
        s_lat = [
            stru.lattice.a,
            stru.lattice.b,
            stru.lattice.c,
            stru.lattice.alpha,
            stru.lattice.beta,
            stru.lattice.gamma,
        ]
        f_lat = [2.489016, 2 * 2.489016, 3 * 2.489016, 60.0, 60.0, 60.0]
        for i in range(len(s_lat)):
            self.assertAlmostEqual(s_lat[i], f_lat[i])
        s_els = [a.element for a in stru]
        self.assertEqual(s_els, 6 * ["Ni"])
        a5 = stru[5]
        s_xyz = a5.xyz
        f_xyz = [0.0, 1.0 / 2.0, 2.0 / 3.0]
        for i in range(3):
            self.assertAlmostEqual(s_xyz[i], f_xyz[i])
        s_o = a5.occupancy
        f_o = 1.0
        self.assertAlmostEqual(s_o, f_o)
        s_U = [a5.U[ij[0], ij[1]] for ij in [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]]
        f_U = 3 * [0.00126651] + 3 * [-0.00042217]
        for i in range(len(s_U)):
            self.assertAlmostEqual(s_U[i], f_U[i])
        return

    def test_read_pdffit_bad(self):
        """check exceptions when reading invalid pdffit file"""
        stru = self.stru
        self.assertRaises(StructureFormatError, stru.read, self.datafile("Ni-bad.stru"), self.format)
        self.assertRaises(StructureFormatError, stru.read, self.datafile("bucky.xyz"), self.format)
        return

    def test_writeStr_pdffit(self):
        """check writing of normal xyz file"""
        stru = self.stru
        stru.read(self.datafile("Ni.stru"), self.format)
        with open(self.datafile("Ni.stru")) as fp:
            f_s = fp.read()
        f_s = re.sub("[ \t]+", " ", f_s)
        f_s = re.sub("[ \t]+\n", "\n", f_s)
        s_s = stru.writeStr(self.format)
        s_s = re.sub("[ \t]+", " ", s_s)
        self.assertEqual(f_s, s_s)
        return

    def test_huge_occupancy(self):
        """check structure with huge occupancy can be read."""
        self.stru.read(self.datafile("Ni.stru"), self.format)
        self.stru[0].occupancy = 16e16
        s_s = self.stru.writeStr(self.format)
        stru1 = Structure()
        stru1.readStr(s_s, self.format)
        self.assertEqual(16e16, stru1[0].occupancy)
        return

    def test_ignored_lines(self):
        """check skipping of ignored lines in the header"""
        r1 = "ignored record 1"
        r2 = "ignored record 2"
        with open(self.datafile("Ni.stru")) as fp:
            ni_lines = fp.readlines()
        ni_lines.insert(2, r1 + "\n")
        ni_lines.insert(4, r2 + "\n")
        s_s1 = "".join(ni_lines)
        p = self.stru.readStr(s_s1, self.format)
        self.assertEqual([r1, r2], p.ignored_lines)
        ni_lines.insert(-3, r1 + "\n")
        s_s2 = "".join(ni_lines)
        self.assertRaises(StructureFormatError, self.stru.readStr, s_s2, self.format)
        return

    def test_spdiameter_parsing(self):
        """check parsing of spdiameter record from a file."""
        stru = self.stru
        stru.read(self.datafile("Ni.stru"), self.format)
        self.assertEqual(0, stru.pdffit["spdiameter"])
        snoshape = stru.writeStr(format=self.format)
        self.assertTrue(not re.search("(?m)^shape", snoshape))
        # produce a string with non-zero spdiameter
        stru.pdffit["spdiameter"] = 13
        s13 = stru.writeStr(format=self.format)
        self.assertTrue(re.search("(?m)^shape +sphere, ", s13))
        stru13 = Structure()
        stru13.readStr(s13)
        self.assertEqual(13, stru13.pdffit["spdiameter"])
        with open(self.datafile("Ni.stru")) as fp:
            ni_lines = fp.readlines()
        ni_lines.insert(3, "shape invalid, 7\n")
        sbad = "".join(ni_lines)
        self.assertRaises(StructureFormatError, self.stru.readStr, sbad, format=self.format)
        return

    def test_stepcut_parsing(self):
        """check parsing of stepcut record from a file."""
        stru = self.stru
        stru.read(self.datafile("Ni.stru"), self.format)
        self.assertEqual(0, stru.pdffit["stepcut"])
        snoshape = stru.writeStr(format=self.format)
        self.assertTrue(not re.search("(?m)^shape", snoshape))
        # produce a string with non-zero stepcut
        stru.pdffit["stepcut"] = 13
        s13 = stru.writeStr(format=self.format)
        self.assertTrue(re.search("(?m)^shape +stepcut, ", s13))
        stru13 = Structure()
        stru13.readStr(s13)
        self.assertEqual(13, stru13.pdffit["stepcut"])
        with open(self.datafile("Ni.stru")) as fp:
            ni_lines = fp.readlines()
        ni_lines.insert(3, "shape invalid, 7\n")
        sbad = "".join(ni_lines)
        self.assertRaises(StructureFormatError, self.stru.readStr, sbad, format=self.format)
        return


# End of class TestP_pdffit

# ----------------------------------------------------------------------------

if __name__ == "__main__":
    unittest.main()

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

"""Unit tests for diffpy.structure.parsers.p_discus module
"""

import re
import unittest

import pytest

from diffpy.structure import Structure
from diffpy.structure.parsers import StructureFormatError

# ----------------------------------------------------------------------------


class TestP_discus(unittest.TestCase):
    """test Parser for PDFFit file format"""

    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    def setUp(self):
        self.stru = Structure()
        self.format = "discus"
        self.places = 8

    def test_read_discus_Ni(self):
        """check reading of Ni structure in discus format"""
        stru = self.stru
        stru.read(self.datafile("Ni-discus.stru"), self.format)
        f_title = "structure Ni  FCC"
        self.assertEqual(f_title, stru.title)
        self.assertEqual("Fm-3m", stru.pdffit["spcgr"])
        # cell record
        abcABG = (3.52, 3.52, 3.52, 90.0, 90.0, 90.0)
        self.assertEqual(abcABG, stru.lattice.abcABG())
        # ncell
        self.assertEqual([1, 1, 1, 4], stru.pdffit["ncell"])
        self.assertEqual(4, len(stru))
        # first atom
        a0 = stru[0]
        self.assertEqual((0.0, 0.0, 0.0), tuple(a0.xyz))
        self.assertTrue(not a0.anisotropy)
        Biso0 = 0.1
        self.assertAlmostEqual(Biso0, a0.Bisoequiv, self.places)
        return

    def test_except_other_formats(self):
        """check exceptions when reading files in other formats"""
        badfiles = [
            "LiCl-bad.cif",
            "PbTe.cif",
            "arginine.pdb",
            "ZnSb_RT_Q28X_VM_20_fxiso.rstr",
            "Ni-bad.stru",
            "Ni.stru",
            "BubbleRaftShort.xcfg",
            "bucky-bad1.xyz",
            "bucky-bad2.xyz",
            "bucky-plain-bad.xyz",
            "bucky-plain.xyz",
            "bucky-raw.xyz",
            "bucky.xyz",
            "hexagon-raw-bad.xyz",
            "hexagon-raw.xyz",
        ]
        for ft in badfiles:
            ff = self.datafile(ft)
            self.assertRaises(StructureFormatError, self.stru.read, ff, format=self.format)
        return

    def test_ignored_lines(self):
        """check skipping of ignored lines in the header"""
        r1 = "ignored record 1\n"
        r2 = "ignored record 2\n"
        with open(self.datafile("Ni-discus.stru")) as fp:
            ni_lines = fp.readlines()
        ni_lines.insert(2, r1)
        ni_lines.insert(4, r2)
        s_s1 = "".join(ni_lines)
        p = self.stru.readStr(s_s1, self.format)
        self.assertEqual([r1.rstrip(), r2.rstrip()], p.ignored_lines)
        ni_lines.append(r1)
        s_s2 = "".join(ni_lines)
        self.assertRaises(StructureFormatError, self.stru.readStr, s_s2, self.format)
        return

    def test_spdiameter_parsing(self):
        """check parsing of spdiameter record from a file."""
        stru = self.stru
        stru.read(self.datafile("Ni-discus.stru"), self.format)
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
        stru.read(self.datafile("Ni-discus.stru"), self.format)
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


# End of class TestP_discus

# ----------------------------------------------------------------------------

if __name__ == "__main__":
    unittest.main()

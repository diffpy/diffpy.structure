#!/usr/bin/env python
##############################################################################
#
# diffpy.structure  by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2006 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas, Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################

"""Unit tests for structure.parsers module.
"""

import os
import re
import tempfile
import unittest

import numpy
import pytest

from diffpy.structure import Atom, Lattice, Structure
from diffpy.structure.structureerrors import StructureFormatError

# ----------------------------------------------------------------------------


class TestP_xyz(unittest.TestCase):
    """test Parser for xyz file format"""

    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    def setUp(self):
        self.stru = Structure()
        self.format = "xyz"
        self.tmpnames = []
        return

    def tearDown(self):
        for f in self.tmpnames:
            os.remove(f)
        return

    def mktmpfile(self):
        with tempfile.NamedTemporaryFile(delete=False) as ftmp:
            self.tmpnames.append(ftmp.name)
        return self.tmpnames[-1]

    def test_read_xyz(self):
        """check reading of normal xyz file"""
        stru = self.stru
        stru.read(self.datafile("bucky.xyz"), self.format)
        s_els = [a.element for a in stru]
        self.assertEqual(stru.title, "bucky-ball")
        self.assertEqual(s_els, 60 * ["C"])
        return

    def test_read_xyz_bad(self):
        """check exceptions when reading invalid xyz file"""
        stru = self.stru
        self.assertRaises(StructureFormatError, stru.read, self.datafile("bucky-bad1.xyz"), self.format)
        self.assertRaises(StructureFormatError, stru.read, self.datafile("bucky-bad2.xyz"), self.format)
        self.assertRaises(StructureFormatError, stru.read, self.datafile("bucky-plain.xyz"), self.format)
        self.assertRaises(StructureFormatError, stru.read, self.datafile("hexagon-raw.xy"), self.format)
        return

    def test_writeStr_xyz(self):
        """check string representation of normal xyz file"""
        stru = self.stru
        stru.title = "test of writeStr"
        stru.lattice = Lattice(1.0, 2.0, 3.0, 90.0, 90.0, 90.0)
        stru[:] = [Atom("H", [1.0, 1.0, 1.0]), Atom("Cl", [3.0, 2.0, 1.0])]
        s1 = stru.writeStr(self.format)
        s1 = re.sub("[ \t]+", " ", s1)
        s0 = "2\n%s\nH 1 2 3\nCl 3 4 3\n" % stru.title
        self.assertEqual(s1, s0)
        return

    def test_write_xyz(self):
        """check writing of normal xyz file"""
        stru = self.stru
        stru.title = "test of writeStr"
        stru.lattice = Lattice(1.0, 2.0, 3.0, 90.0, 90.0, 90.0)
        stru[:] = [Atom("H", [1.0, 1.0, 1.0]), Atom("Cl", [3.0, 2.0, 1.0])]
        filename = self.mktmpfile()
        stru.write(filename, self.format)
        with open(filename) as fp:
            f_s = fp.read()
        f_s = re.sub("[ \t]+", " ", f_s)
        s_s = "2\n%s\nH 1 2 3\nCl 3 4 3\n" % stru.title
        self.assertEqual(f_s, s_s)
        return


# End of class TestP_xyz

# ----------------------------------------------------------------------------


class TestP_rawxyz(unittest.TestCase):
    """test Parser for rawxyz file format"""

    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    def setUp(self):
        self.stru = Structure()
        self.format = "rawxyz"
        return

    def test_read_plainxyz(self):
        """check reading of a plain xyz file"""
        stru = self.stru
        stru.read(self.datafile("bucky-plain.xyz"), self.format)
        s_els = [a.element for a in stru]
        self.assertEqual(stru.title, "bucky-plain")
        self.assertEqual(s_els, 60 * ["C"])
        return

    def test_read_plainxyz_bad(self):
        """check exceptions when reading invalid plain xyz file"""
        stru = self.stru
        self.assertRaises(StructureFormatError, stru.read, self.datafile("bucky-plain-bad.xyz"), self.format)
        return

    def test_read_rawxyz(self):
        """check reading of raw xyz file"""
        stru = self.stru
        stru.read(self.datafile("bucky-raw.xyz"), self.format)
        s_els = [a.element for a in stru]
        self.assertEqual(stru.title, "bucky-raw")
        self.assertEqual(s_els, 60 * [""])
        stru.read(self.datafile("hexagon-raw.xyz"), self.format)
        zs = [a.xyz[-1] for a in stru]
        self.assertEqual(zs, 6 * [0.0])
        return

    def test_read_rawxyz_bad(self):
        """check exceptions when reading unsupported xy file"""
        stru = self.stru
        self.assertRaises(StructureFormatError, stru.read, self.datafile("hexagon-raw-bad.xyz"), self.format)
        self.assertRaises(StructureFormatError, stru.read, self.datafile("hexagon-raw.xy"), self.format)
        return

    def test_writeStr_rawxyz(self):
        """check writing of normal xyz file"""
        stru = self.stru
        stru.title = "test of writeStr"
        stru.lattice = Lattice(1.0, 2.0, 3.0, 90.0, 90.0, 90.0)
        # plain version
        stru[:] = [Atom("H", [1.0, 1.0, 1.0])]
        s1 = stru.writeStr(self.format)
        s1 = re.sub("[ \t]+", " ", s1)
        s0 = "H 1 2 3\n"
        # brutal raw version
        stru[0].element = ""
        s1 = stru.writeStr(self.format)
        s0 = "1 2 3\n"
        self.assertEqual(s1, s0)
        return


# End of class TestP_rawxyz

# ----------------------------------------------------------------------------


class TestP_pdb(unittest.TestCase):
    """test Parser for PDB file format"""

    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    def setUp(self):
        self.stru = Structure()
        self.format = "pdb"
        self.places = 3

    def test_read_pdb_arginine(self):
        """check reading of arginine PDB file"""
        stru = self.stru
        stru.read(self.datafile("arginine.pdb"), self.format)
        f_els = [
            "N",
            "C",
            "C",
            "O",
            "C",
            "C",
            "C",
            "N",
            "C",
            "N",
            "N",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
            "H",
            "O",
            "H",
        ]
        s_els = [a.element for a in stru]
        self.assertEqual(s_els, f_els)
        s_lat = [
            stru.lattice.a,
            stru.lattice.b,
            stru.lattice.c,
            stru.lattice.alpha,
            stru.lattice.beta,
            stru.lattice.gamma,
        ]
        f_lat = [1.0, 1.0, 1.0, 90.0, 90.0, 90.0]
        self.assertEqual(s_lat, f_lat)
        a0 = stru[0]
        self.assertTrue(numpy.allclose(a0.xyz, [0.735, 2.219, 1.389]))

    def test_rwStr_pdb_CdSe(self):
        """check conversion to PDB file format"""
        stru = self.stru
        stru.read(self.datafile("CdSe_bulk.stru"), "pdffit")
        s = stru.writeStr(self.format)
        # all lines should be 80 characters long
        linelens = [len(line) for line in s.split("\n") if line != ""]
        self.assertEqual(linelens, len(linelens) * [80])
        # now clean and re-read structure
        stru = Structure()
        stru.readStr(s, self.format)
        s_els = [a.element for a in stru]
        f_els = ["Cd", "Cd", "Se", "Se"]
        self.assertEqual(s_els, f_els)
        s_lat = [
            stru.lattice.a,
            stru.lattice.b,
            stru.lattice.c,
            stru.lattice.alpha,
            stru.lattice.beta,
            stru.lattice.gamma,
        ]
        f_lat = [4.235204, 4.235204, 6.906027, 90.0, 90.0, 120.0]
        self.assertTrue(numpy.allclose(s_lat, f_lat, atol=5e-4))
        a0 = stru[0]
        s_Uii = [a0.U[i, i] for i in range(3)]
        f_Uii = [0.01303035, 0.01303035, 0.01401959]
        self.assertTrue(numpy.allclose(s_Uii, f_Uii, atol=5e-4))
        s_sigUii = [a0.sigU[i, i] for i in range(3)]
        f_sigUii = [0.00011127, 0.00011127, 0.00019575]
        self.assertTrue(numpy.allclose(s_sigUii, f_sigUii, atol=5e-4))
        s_title = stru.title
        f_title = "Cell structure file of CdSe #186"
        self.assertEqual(s_title, f_title)


# End of class TestP_pdb

# ----------------------------------------------------------------------------


class TestP_xcfg(unittest.TestCase):
    """test Parser for XCFG file format"""

    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    def setUp(self):
        self.stru = Structure()
        self.format = "xcfg"
        self.places = 6

    def test_read_xcfg(self):
        """check reading of BubbleRaft XCFG file"""
        stru = self.stru
        stru.read(self.datafile("BubbleRaftShort.xcfg"), self.format)
        f_els = 500 * ["Ar"]
        s_els = [a.element for a in stru]
        self.assertEqual(s_els, f_els)
        self.assertAlmostEqual(stru.distance(82, 357), 47.5627, 3)
        s_lat = [
            stru.lattice.a,
            stru.lattice.b,
            stru.lattice.c,
            stru.lattice.alpha,
            stru.lattice.beta,
            stru.lattice.gamma,
        ]
        f_lat = [127.5, 119.5, 3.0, 90.0, 90.0, 90.0]
        self.assertTrue(numpy.allclose(s_lat, f_lat))
        return

    def test_rwStr_xcfg_CdSe(self):
        """check conversion to XCFG file format"""
        stru = self.stru
        stru.read(self.datafile("CdSe_bulk.stru"), "pdffit")
        s = stru.writeStr(self.format)
        stru = Structure()
        stru.readStr(s, self.format)
        s_els = [a.element for a in stru]
        f_els = ["Cd", "Cd", "Se", "Se"]
        self.assertEqual(s_els, f_els)
        s_lat = [
            stru.lattice.a,
            stru.lattice.b,
            stru.lattice.c,
            stru.lattice.alpha,
            stru.lattice.beta,
            stru.lattice.gamma,
        ]
        f_lat = [4.235204, 4.235204, 6.906027, 90.0, 90.0, 120.0]
        self.assertTrue(numpy.allclose(s_lat, f_lat))
        a0 = stru[0]
        s_Uii = [a0.U[i, i] for i in range(3)]
        f_Uii = [0.01303035, 0.01303035, 0.01401959]
        self.assertTrue(numpy.allclose(s_Uii, f_Uii))


# End of class TestP_xcfg

# ----------------------------------------------------------------------------

if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python

"""Unit tests for the loadStructure factory.
"""

import unittest

import pytest

from diffpy.structure import PDFFitStructure, Structure, loadStructure
from diffpy.structure.structureerrors import StructureFormatError


##############################################################################
class TestLoadStructure(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.datafile = datafile

    def test_xcfg(self):
        """check loading of atomeye xcfg format"""
        f = self.datafile("BubbleRaftShort.xcfg")
        stru = loadStructure(f)
        self.assertTrue(type(stru) is Structure)
        self.assertRaises(StructureFormatError, loadStructure, f, "xyz")
        return

    def test_discus(self):
        """check loading of discus file format"""
        f = self.datafile("Ni-discus.stru")
        stru = loadStructure(f)
        self.assertTrue(type(stru) is PDFFitStructure)
        return

    def test_cif(self):
        """check loading of CIF file format"""
        f = self.datafile("PbTe.cif")
        stru = loadStructure(f)
        self.assertTrue(isinstance(stru, Structure))
        self.assertFalse(isinstance(stru, PDFFitStructure))
        return

    def test_badfile(self):
        """check loading of CIF file format"""
        f = self.datafile("Ni-bad.stru")
        self.assertRaises(StructureFormatError, loadStructure, f)
        return

    def test_goodkwarg(self):
        """check loading of CIF file and passing of parser keyword argument."""
        f = self.datafile("graphite.cif")
        stru = loadStructure(f, eps=1e-10)
        self.assertEqual(8, len(stru))
        return

    def test_badkwarg(self):
        """check loading of xyz file format with invalid keyword argument"""
        f = self.datafile("bucky.xyz")
        self.assertRaises(TypeError, loadStructure, f, eps=1e-10)
        return


# End of class TestLoadStructure

# ----------------------------------------------------------------------------

if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python

"""Unit tests for the loadStructure factory.
"""

import unittest
from diffpy.Structure.tests.testutils import datafile
from diffpy.Structure import loadStructure
from diffpy.Structure import Structure, PDFFitStructure, StructureFormatError


##############################################################################
class TestLoadStructure(unittest.TestCase):

    def test_xcfg(self):
        """check loading of atomeye xcfg format
        """
        f = datafile('BubbleRaftShort.xcfg')
        stru = loadStructure(f)
        self.failUnless(type(stru) is Structure)
        self.assertRaises(StructureFormatError,
                loadStructure, f, 'xyz')
        return

    def test_discus(self):
        """check loading of discus file format
        """
        f = datafile('Ni-discus.stru')
        stru = loadStructure(f)
        self.failUnless(type(stru) is PDFFitStructure)
        return

    def test_cif(self):
        """check loading of CIF file format
        """
        f = datafile('PbTe.cif')
        stru = loadStructure(f)
        self.failUnless(isinstance(stru, Structure))
        self.failIf(isinstance(stru, PDFFitStructure))
        return

    def test_badfile(self):
        """check loading of CIF file format
        """
        f = datafile('Ni-bad.stru')
        self.assertRaises(StructureFormatError, loadStructure, f)
        return


# End of class TestLoadStructure

if __name__ == '__main__':
    unittest.main()

# End of file

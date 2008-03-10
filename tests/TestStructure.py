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

"""Unit tests for Structure class.
"""

__id__ = "$Id$"

import os
import unittest

# useful variables
thisfile = locals().get('__file__', 'TestStructure.py')
tests_dir = os.path.dirname(os.path.abspath(thisfile))
testdata_dir = os.path.join(tests_dir, 'testdata')

from diffpy.Structure import Structure, StructureFormatError
from diffpy.Structure import Lattice
from diffpy.Structure import Atom

##############################################################################
class TestStructure(unittest.TestCase):
    """test methods of Structure class"""

    def setUp(self):
        self.stru = Structure( [ Atom('C', [0,0,0]), Atom('C', [1,1,1]) ],
                lattice=Lattice(1, 1, 1, 90, 90, 120) )
        self.places = 12

    def assertListAlmostEqual(self, l1, l2, places=None):
        """wrapper for list comparison"""
        if places is None: places = self.places
        self.assertEqual(len(l1), len(l2))
        for i in range(len(l1)):
            self.assertAlmostEqual(l1[i], l2[i], places)

    # FIXME move into TestAtom
    def test_cartesian(self):
        """check conversion to cartesian coordinates"""
        from math import sqrt
        stru = self.stru
        s_rc0 = stru[0].xyz_cartn
        f_rc0 = 3*[0.0]
        s_rc1 = stru[1].xyz_cartn
        f_rc1 = [sqrt(0.75), 0.5, 1.]
        self.assertListAlmostEqual(s_rc0, f_rc0)
        self.assertListAlmostEqual(s_rc1, f_rc1)

#   def test___init__(self):
#       """check Structure.__init__()
#       """
#       return
#
#   def test___str__(self):
#       """check Structure.__str__()
#       """
#       return
#
#   def test_addNewAtom(self):
#       """check Structure.addNewAtom()
#       """
#       return
#
#   def test_getLastAtom(self):
#       """check Structure.getLastAtom()
#       """
#       return

    def test_getAtom(self):
        """check Structure.getAtom()
        """
        a0 = self.stru[0]
        a1 = self.stru[1]
        # check execeptions for invalid arguments
        self.assertRaises(ValueError, self.stru.getAtom, 300)
        self.assertRaises(ValueError, self.stru.getAtom, -44)
        self.assertRaises(ValueError, self.stru.getAtom, "Na")
        # check returned values
        self.failUnless(a0 is self.stru.getAtom(0))
        self.failUnless(a1 is self.stru.getAtom(1))
        self.failUnless(a0 is self.stru.getAtom("C1"))
        self.failUnless(a1 is self.stru.getAtom("C2"))
        # check if labels get properly updated
        cdsefile = os.path.join(testdata_dir, 'CdSe_bulk.stru')
        cdse = Structure(filename=cdsefile)
        self.stru[1:1] = cdse
        self.failUnless(a0 is self.stru.getAtom("C1"))
        self.failUnless(a1 is self.stru.getAtom("C2"))
        self.failUnless(self.stru[1] is self.stru.getAtom("Cd1"))
        return


    def test_getLabels(self):
        """check Structure.getLabels()
        """
        self.assertEqual(["C1", "C2"], self.stru.getLabels())
        pbtefile = os.path.join(testdata_dir, 'PbTe.cif')
        self.stru.read(pbtefile, format='cif')
        labels = self.stru.getLabels()
        self.assertEqual("Pb2+1", labels[0])
        self.assertEqual("Pb2+4", labels[3])
        self.assertEqual("Te1", labels[4])
        self.assertEqual("Te4", labels[-1])
        return


    def test_distance(self):
        """check Structure.distance()
        """
        from math import sqrt
        self.assertRaises(ValueError, self.stru.distance, 333, "C1")
        self.assertRaises(ValueError, self.stru.distance, "C", "C1")
        self.assertAlmostEqual(sqrt(2.0),
                self.stru.distance(0, 1), self.places)
        self.assertAlmostEqual(sqrt(2.0),
                self.stru.distance("C1", "C2"), self.places)
        self.assertEqual(0, self.stru.distance(0, "C1"))
        return

#   def test_angle(self):
#       """check Structure.angle()
#       """
#       return

    def test_placeInLattice(self):
        """check Structure.placeInLattice() -- conversion of coordinates
        """
        stru = self.stru
        new_lattice = Lattice(.5, .5, .5, 90, 90, 60)
        stru.placeInLattice(new_lattice)
        a0 = stru[0]
        self.assertListAlmostEqual(a0.xyz, 3*[0.0])
        a1 = stru[1]
        self.assertListAlmostEqual(a1.xyz, [2.0, 0.0, 2.0])

#   def test_read(self):
#       """check Structure.read()
#       """
#       return
#
#   def test_readStr(self):
#       """check Structure.readStr()
#       """
#       return
#
#   def test_write(self):
#       """check Structure.write()
#       """
#       return
#
#   def test_writeStr(self):
#       """check Structure.writeStr()
#       """
#       return
#
#   def test_append(self):
#       """check Structure.append()
#       """
#       return
#
#   def test_insert(self):
#       """check Structure.insert()
#       """
#       return
#
#   def test_extend(self):
#       """check Structure.extend()
#       """
#       return
#
#   def test___setitem__(self):
#       """check Structure.__setitem__()
#       """
#       return
#
#   def test___setslice__(self):
#       """check Structure.__setslice__()
#       """
#       return
#
#   def test__get_lattice(self):
#       """check Structure._get_lattice()
#       """
#       return
#
#   def test__set_lattice(self):
#       """check Structure._set_lattice()
#       """
#       return
#
#   def test__update_labels(self):
#       """check Structure._update_labels()
#       """
#       return
#
#   def test__uncache(self):
#       """check Structure._uncache()
#       """
#       return

# End of class TestStructure

if __name__ == '__main__':
    unittest.main()

# End of file

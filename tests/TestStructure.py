##############################################################################
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
##############################################################################

"""Unit tests for Structure class.
"""

__id__ = "$Id$"

import os
import copy
import unittest
import numpy

# useful variables
thisfile = locals().get('__file__', 'TestStructure.py')
tests_dir = os.path.dirname(os.path.abspath(thisfile))
testdata_dir = os.path.join(tests_dir, 'testdata')
cdsefile = os.path.join(testdata_dir, 'CdSe_bulk.stru')
pbtefile = os.path.join(testdata_dir, 'PbTe.cif')

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

    def test___copy__(self):
        """check Structure.__copy__()
        """
        cdse = Structure(filename=cdsefile)
        cdse_str = cdse.writeStr('pdffit')
        cdse2 = copy.copy(cdse)
        self.assertEqual(cdse_str, cdse2.writeStr('pdffit'))
        self.failIf(cdse.lattice is cdse2.lattice)
        sameatoms = set(cdse).intersection(cdse2)
        self.failIf(sameatoms)
        return

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


    def test_angle(self):
        """check Structure.angle()
        """
        cdse = Structure(filename=cdsefile)
        self.assertEqual(109, round(cdse.angle(0, 2, 1)))
        self.assertEqual(109, round(cdse.angle("Cd1", "Se1", "Cd2")))
        return


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

    def test_aslist(self):
        """check Structure.aslist()
        """
        lst = self.stru.aslist()
        self.assertEqual(tuple(lst), tuple(self.stru))
        self.assertEqual(list, type(lst))
        return


    def test_append(self):
        """check Structure.append()
        """
        a = Atom("Si", (0.1, 0.2, 0.3))
        lat = self.stru.lattice
        self.stru.append(a)
        alast = self.stru[-1]
        self.assertEqual(3, len(self.stru))
        self.assertEqual('Si', alast.element)
        self.failUnless(lat is alast.lattice)
        self.failUnless(numpy.array_equal(a.xyz, alast.xyz))
        self.failIf(a is alast)
        self.failIf(lat is a.lattice)
        return


    def test_insert(self):
        """check Structure.insert()
        """
        a = Atom("Si", (0.1, 0.2, 0.3))
        lat = self.stru.lattice
        self.stru.insert(1, a)
        a1 = self.stru[1]
        self.assertEqual(3, len(self.stru))
        self.assertEqual('Si', a1.element)
        self.failUnless(lat is a1.lattice)
        self.failUnless(numpy.array_equal(a.xyz, a1.xyz))
        self.failIf(a is a1)
        self.failIf(lat is a.lattice)
        return


    def test_extend(self):
        """check Structure.extend()
        """
        stru = self.stru
        cdse = Structure(filename=cdsefile)
        lst = stru.aslist()
        stru.extend(cdse)
        self.assertEqual(6, len(stru))
        self.failUnless(all([a.lattice is stru.lattice for a in stru]))
        self.failUnless(stru.lattice is a.lattice)
        self.assertEqual(lst, stru.aslist()[:2])
        self.assertNotEqual(stru[-1], cdse[-1])
        return


    def test___getitem__(self):
        """check Structure.__getitem__()
        """
        stru = self.stru
        self.failUnless(stru[0] is stru.aslist()[0])
        intidx = range(len(stru))[::-1]
        self.assertEqual(stru[intidx].aslist(), stru.aslist()[::-1])
        flagidx = (numpy.arange(len(stru)) > 0)
        self.assertEqual(stru[flagidx].aslist(), stru.aslist()[1:])
        return


    def test___setitem__(self):
        """check Structure.__setitem__()
        """
        a = Atom("Si", (0.1, 0.2, 0.3))
        lat = self.stru.lattice
        self.stru[1] = a
        a1 = self.stru[1]
        self.assertEqual(2, len(self.stru))
        self.assertEqual('Si', a1.element)
        self.failUnless(lat is a1.lattice)
        self.failUnless(numpy.array_equal(a.xyz, a1.xyz))
        self.failIf(a is a1)
        self.failIf(lat is a.lattice)
        return


    def test___setslice__(self):
        """check Structure.__setslice__()
        """
        a = Atom("Si", (0.1, 0.2, 0.3))
        lat = self.stru.lattice
        self.stru[:] = [a]
        a0 = self.stru[0]
        self.assertEqual(1, len(self.stru))
        self.assertEqual('Si', a0.element)
        self.failUnless(lat is a0.lattice)
        self.failUnless(numpy.array_equal(a.xyz, a0.xyz))
        self.failIf(a is a0)
        self.failIf(lat is a.lattice)
        return


    def test___add__(self):
        """check Structure.__add__()
        """
        stru = self.stru
        cdse = Structure(filename=cdsefile)
        total = stru + cdse
        self.assertEqual(6, len(total))
        ta0 = total[0]
        tam1 = total[-1]
        self.assertEqual('C', ta0.element)
        self.failUnless(numpy.array_equal(stru[0].xyz, ta0.xyz))
        self.assertEqual('Se', tam1.element)
        self.failUnless(numpy.array_equal(cdse[-1].xyz, tam1.xyz))
        self.failIf(total.lattice in (stru.lattice, cdse.lattice))
        self.failUnless(all([a.lattice is total.lattice for a in total]))
        return


    def test___iadd__(self):
        """check Structure.__iadd__()
        """
        stru = self.stru
        lat0 = stru.lattice
        lst = stru.aslist()
        cdse = Structure(filename=cdsefile)
        stru += cdse
        self.assertEqual(6, len(stru))
        self.assertEqual(lst, stru[:2].aslist())
        am1 = stru[-1]
        self.assertEqual('Se', am1.element)
        self.failUnless(numpy.array_equal(cdse[-1].xyz, am1.xyz))
        self.failUnless(lat0 is stru.lattice)
        self.failIf(stru.lattice is cdse.lattice)
        self.failUnless(all([a.lattice is stru.lattice for a in stru]))
        return


    def test___sub__(self):
        """check Structure.__sub__()
        """
        cdse = Structure(filename=cdsefile)
        cadmiums = cdse - cdse[2:]
        self.assertEqual(2, len(cadmiums))
        self.assertEqual('Cd', cadmiums[0].element)
        self.assertEqual('Cd', cadmiums[1].element)
        self.failUnless(numpy.array_equal(cdse[0].xyz, cadmiums[0].xyz))
        self.failUnless(numpy.array_equal(cdse[1].xyz, cadmiums[1].xyz))
        self.failIf(cdse[0] is cadmiums[0])
        self.failIf(cdse.lattice is cadmiums.lattice)
        return


    def test___isub__(self):
        """check Structure.__isub__()
        """
        cdse = Structure(filename=cdsefile)
        lat = cdse.lattice
        lst = cdse.aslist()
        cdse -= cdse[2:]
        self.assertEqual(2, len(cdse))
        self.assertEqual(4, len(lst))
        self.assertEqual('Cd', cdse[0].element)
        self.assertEqual('Cd', cdse[1].element)
        self.assertEqual(lat, cdse.lattice)
        self.assertEqual(lst[:2], cdse.aslist())
        return


    def test___mul__(self):
        """check Structure.__mul__()
        """
        cdse = Structure(filename=cdsefile)
        self.assertEqual(12, len(set(3 * cdse)))
        self.assertEqual(12, len(set(cdse * 3)))
        cdsex3 = 3 * cdse
        self.assertEqual(12, len(cdsex3))
        self.assertEqual(3 * 'Cd Cd Se Se'.split(),
            [a.element for a in cdsex3])
        self.failUnless(numpy.array_equal(3 * [a.xyz for a in cdse],
            [a.xyz for a in cdsex3]))
        self.failIf(set(cdse).intersection(cdsex3))
        self.failIf(cdse.lattice is cdsex3.lattice)
        return


    def test___imul__(self):
        """check Structure.__imul__()
        """
        cdse = Structure(filename=cdsefile)
        lat = cdse.lattice
        els = cdse.element
        xyz = cdse.xyz
        lst = cdse.aslist()
        cdse *= 2
        self.assertEqual(8, len(cdse))
        self.assertEqual(lst, cdse[:4].aslist())
        self.assertEqual(numpy.tile(els, 2).tolist(), cdse.element.tolist())
        self.failUnless(numpy.array_equal(numpy.tile(xyz, (2, 1)), cdse.xyz))
        self.assertEqual(8, len(set(cdse)))
        self.assertEqual(8 * [lat], [a.lattice for a in cdse])
        return

    def test__get_lattice(self):
        """check Structure._get_lattice()
        """
        lat = Lattice()
        stru = Structure()
        self.assertEqual((1, 1, 1, 90, 90, 90), stru.lattice.abcABG())
        stru2 = Structure(lattice=lat)
        self.failUnless(lat is stru2.lattice)
        return

    def test__set_lattice(self):
        """check Structure._set_lattice()
        """
        lat = Lattice()
        self.stru.lattice = lat
        self.assertEqual(2 * [lat], [a.lattice for a in self.stru])
        return

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

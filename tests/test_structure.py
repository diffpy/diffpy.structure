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

"""Unit tests for Structure class.
"""


import copy
import pickle
import unittest

import numpy
import pytest

from diffpy.structure import Atom, Lattice, Structure

# ----------------------------------------------------------------------------


class TestStructure(unittest.TestCase):
    """test methods of Structure class"""

    @pytest.fixture(autouse=True)
    def prepare_fixture(self, datafile):
        self.cdsefile = datafile("CdSe_bulk.stru")
        self.teifile = datafile("TeI.cif")
        self.pbtefile = datafile("PbTe.cif")

    _loaded_structures = {}

    def setUp(self):
        self.stru = Structure([Atom("C", [0, 0, 0]), Atom("C", [1, 1, 1])], lattice=Lattice(1, 1, 1, 90, 90, 120))
        # useful variables

        if not self._loaded_structures:
            self._loaded_structures.update(
                [
                    ("cdse", Structure(filename=self.cdsefile)),
                    ("tei", Structure(filename=self.teifile)),
                    ("pbte", Structure(filename=self.pbtefile)),
                ]
            )
        self.__dict__.update(self._loaded_structures)
        self.places = 12
        return

    def test___init__(self):
        """check Structure.__init__()"""
        atoms = [Atom("C", [0, 0, 0]), Atom("C", [0.5, 0.5, 0.5])]
        self.assertRaises(ValueError, Structure, atoms, filename=self.teifile)
        self.assertRaises(ValueError, Structure, lattice=Lattice(), filename=self.teifile)
        self.assertRaises(ValueError, Structure, title="test", filename=self.teifile)
        stru1 = Structure(title="my title")
        self.assertEqual("my title", stru1.title)
        stru2a = Structure(atoms)
        stru2b = Structure(iter(atoms))
        stru2c = Structure(a for a in atoms)
        s2a = str(stru2a)
        self.assertEqual(s2a, str(stru2b))
        self.assertEqual(s2a, str(stru2c))
        return

    def test_copy(self):
        """check Structure.copy()"""

        class MyDerivedStructure(Structure):
            def __copy__(self):
                rv = MyDerivedStructure()
                Structure.__copy__(self, target=rv)
                return rv

            pass

        pbte = self.pbte
        pbte2 = pbte.copy()
        self.assertFalse(pbte2.lattice is pbte.lattice)
        self.assertTrue(numpy.array_equal(pbte.xyz_cartn, pbte2.xyz_cartn))
        self.assertTrue(numpy.array_equal(pbte.U, pbte2.U))
        stru = MyDerivedStructure()
        stru += pbte2[pbte2.element.startswith("Pb")]
        pb3 = stru.copy()
        self.assertTrue(isinstance(pb3, MyDerivedStructure))
        self.assertTrue(all(pb3.element == "Pb2+"))
        self.assertEqual(4, len(pb3))
        return

    def test___copy__(self):
        """check Structure.__copy__()"""
        cdse = Structure(filename=self.cdsefile)
        cdse_str = cdse.writeStr("pdffit")
        cdse2 = copy.copy(cdse)
        self.assertEqual(cdse_str, cdse2.writeStr("pdffit"))
        self.assertFalse(cdse.lattice is cdse2.lattice)
        sameatoms = set(cdse).intersection(cdse2)
        self.assertFalse(sameatoms)
        return

    # def test___str__(self):
    #     """check Structure.__str__()"""
    #     return

    # def test_addNewAtom(self):
    #     """check Structure.addNewAtom()"""
    #     return

    # def test_getLastAtom(self):
    #     """check Structure.getLastAtom()"""
    #     return

    def test_assignUniqueLabels(self):
        """check Structure.assignUniqueLabels()"""
        self.assertEqual("", "".join([a.label for a in self.stru]))
        self.stru.assignUniqueLabels()
        self.assertEqual("C1", self.stru[0].label)
        self.assertEqual("C2", self.stru[1].label)
        return

    def test_distance(self):
        """check Structure.distance()"""
        from math import sqrt

        self.stru.assignUniqueLabels()
        self.assertRaises(IndexError, self.stru.distance, 333, "C1")
        self.assertRaises(IndexError, self.stru.distance, "C", "C1")
        self.assertAlmostEqual(sqrt(2.0), self.stru.distance(0, 1), self.places)
        self.assertAlmostEqual(sqrt(2.0), self.stru.distance("C1", "C2"), self.places)
        self.assertEqual(0, self.stru.distance(0, "C1"))
        return

    def test_angle(self):
        """check Structure.angle()"""
        cdse = Structure(filename=self.cdsefile)
        cdse.assignUniqueLabels()
        self.assertEqual(109, round(cdse.angle(0, 2, 1)))
        self.assertEqual(109, round(cdse.angle("Cd1", "Se1", "Cd2")))
        return

    def test_placeInLattice(self):
        """check Structure.placeInLattice() -- conversion of coordinates"""
        stru = self.stru
        new_lattice = Lattice(0.5, 0.5, 0.5, 90, 90, 60)
        stru.placeInLattice(new_lattice)
        a0 = stru[0]
        self.assertTrue(numpy.allclose(a0.xyz, [0.0, 0.0, 0.0]))
        a1 = stru[1]
        self.assertTrue(numpy.allclose(a1.xyz, [2.0, 0.0, 2.0]))

    # def test_read(self):
    #     """check Structure.read()"""
    #     return

    # def test_readStr(self):
    #     """check Structure.readStr()"""
    #     return

    # def test_write(self):
    #     """check Structure.write()"""
    #     return

    # def test_writeStr(self):
    #     """check Structure.writeStr()"""
    #     return

    def test_aslist(self):
        """check Structure.tolist()"""
        lst = self.stru.tolist()
        self.assertEqual(tuple(lst), tuple(self.stru))
        self.assertEqual(list, type(lst))
        return

    def test_append(self):
        """check Structure.append()"""
        a = Atom("Si", (0.1, 0.2, 0.3))
        lat = self.stru.lattice
        self.stru.append(a)
        alast = self.stru[-1]
        self.assertEqual(3, len(self.stru))
        self.assertEqual("Si", alast.element)
        self.assertTrue(lat is alast.lattice)
        self.assertTrue(numpy.array_equal(a.xyz, alast.xyz))
        self.assertFalse(a is alast)
        self.assertFalse(lat is a.lattice)
        return

    def test_insert(self):
        """check Structure.insert()"""
        a = Atom("Si", (0.1, 0.2, 0.3))
        lat = self.stru.lattice
        self.stru.insert(1, a)
        a1 = self.stru[1]
        self.assertEqual(3, len(self.stru))
        self.assertEqual("Si", a1.element)
        self.assertTrue(lat is a1.lattice)
        self.assertTrue(numpy.array_equal(a.xyz, a1.xyz))
        self.assertFalse(a is a1)
        self.assertFalse(lat is a.lattice)
        return

    def test_extend(self):
        """check Structure.extend()"""
        stru = self.stru
        cdse = Structure(filename=self.cdsefile)
        lst = stru.tolist()
        stru.extend(cdse)
        self.assertEqual(6, len(stru))
        self.assertTrue(all(a.lattice is stru.lattice for a in stru))
        self.assertEqual(lst, stru.tolist()[:2])
        self.assertFalse(stru[-1] is cdse[-1])
        return

    def test___getitem__(self):
        """check Structure.__getitem__()"""
        stru = self.stru
        self.assertTrue(stru[0] is stru.tolist()[0])
        intidx = list(range(len(stru)))[::-1]
        self.assertEqual(stru[intidx].tolist(), stru.tolist()[::-1])
        flagidx = numpy.arange(len(stru)) > 0
        self.assertEqual(stru[flagidx].tolist(), stru.tolist()[1:])
        cdse = Structure(self.cdse)
        self.assertEqual([cdse[0], cdse[-2]], cdse[0, -2].tolist())
        cdse013 = cdse.tolist()
        cdse013.pop(2)
        self.assertEqual(cdse013, cdse[:2, 3].tolist())
        self.assertRaises(IndexError, cdse.__getitem__, "Cd1")
        cdse.assignUniqueLabels()
        self.assertTrue(cdse[0] is cdse["Cd1"])
        cdse[0].label = "Hohenzollern"
        self.assertRaises(IndexError, cdse.__getitem__, "Cd1")
        self.assertTrue(cdse[0] is cdse["Hohenzollern"])
        self.assertEqual([cdse[0], cdse[3], cdse[1]], cdse["Hohenzollern", 3:0:-2].tolist())
        stru.label = ["A", "B"]
        self.assertTrue(stru[0] is stru["A"])
        self.assertTrue(stru[1] is stru["B"])
        stru[1].label = "A"
        self.assertRaises(IndexError, stru.__getitem__, "A")
        return

    def test___getitem__slice(self):
        """check Structure.__getitem__() with a slice argument"""
        stru = self.stru
        self.assertEqual([stru[0]], stru[:1].tolist())
        self.assertEqual([stru[1], stru[0]], stru[::-1].tolist())
        return

    def test___setitem__(self):
        """check Structure.__setitem__()"""
        a = Atom("Si", (0.1, 0.2, 0.3))
        lat = self.stru.lattice
        self.stru[1] = a
        a1 = self.stru[1]
        self.assertEqual(2, len(self.stru))
        self.assertEqual("Si", a1.element)
        self.assertTrue(lat is a1.lattice)
        self.assertTrue(numpy.array_equal(a.xyz, a1.xyz))
        self.assertFalse(a is a1)
        self.assertFalse(lat is a.lattice)
        return

    def test___setitem__slice(self):
        """check Structure.__setitem__() with a slice argument"""
        a = Atom("Si", (0.1, 0.2, 0.3))
        lat = self.stru.lattice
        self.stru[:] = [a]
        a0 = self.stru[0]
        self.assertEqual(1, len(self.stru))
        self.assertEqual("Si", a0.element)
        self.assertTrue(lat is a0.lattice)
        self.assertTrue(numpy.array_equal(a.xyz, a0.xyz))
        self.assertFalse(a is a0)
        self.assertFalse(lat is a.lattice)
        return

    def test___add__(self):
        """check Structure.__add__()"""
        stru = self.stru
        cdse = Structure(filename=self.cdsefile)
        total = stru + cdse
        self.assertEqual(6, len(total))
        ta0 = total[0]
        tam1 = total[-1]
        self.assertEqual("C", ta0.element)
        self.assertTrue(numpy.array_equal(stru[0].xyz, ta0.xyz))
        self.assertEqual("Se", tam1.element)
        self.assertTrue(numpy.array_equal(cdse[-1].xyz, tam1.xyz))
        self.assertFalse(total.lattice in (stru.lattice, cdse.lattice))
        self.assertTrue(all([a.lattice is total.lattice for a in total]))
        return

    def test___iadd__(self):
        """check Structure.__iadd__()"""
        stru = self.stru
        lat0 = stru.lattice
        lst = stru.tolist()
        cdse = Structure(filename=self.cdsefile)
        stru += cdse
        self.assertEqual(6, len(stru))
        self.assertEqual(lst, stru[:2].tolist())
        am1 = stru[-1]
        self.assertEqual("Se", am1.element)
        self.assertTrue(numpy.array_equal(cdse[-1].xyz, am1.xyz))
        self.assertTrue(lat0 is stru.lattice)
        self.assertFalse(stru.lattice is cdse.lattice)
        self.assertTrue(all([a.lattice is stru.lattice for a in stru]))
        return

    def test___sub__(self):
        """check Structure.__sub__()"""
        cdse = Structure(filename=self.cdsefile)
        cadmiums = cdse - cdse[2:]
        self.assertEqual(2, len(cadmiums))
        self.assertEqual("Cd", cadmiums[0].element)
        self.assertEqual("Cd", cadmiums[1].element)
        self.assertTrue(numpy.array_equal(cdse[0].xyz, cadmiums[0].xyz))
        self.assertTrue(numpy.array_equal(cdse[1].xyz, cadmiums[1].xyz))
        self.assertFalse(cdse[0] is cadmiums[0])
        self.assertFalse(cdse.lattice is cadmiums.lattice)
        return

    def test___isub__(self):
        """check Structure.__isub__()"""
        cdse = Structure(filename=self.cdsefile)
        lat = cdse.lattice
        lst = cdse.tolist()
        cdse -= cdse[2:]
        self.assertEqual(2, len(cdse))
        self.assertEqual(4, len(lst))
        self.assertEqual("Cd", cdse[0].element)
        self.assertEqual("Cd", cdse[1].element)
        self.assertEqual(lat, cdse.lattice)
        self.assertEqual(lst[:2], cdse.tolist())
        return

    def test___mul__(self):
        """check Structure.__mul__()"""
        cdse = Structure(filename=self.cdsefile)
        self.assertEqual(12, len(set(3 * cdse)))
        self.assertEqual(12, len(set(cdse * 3)))
        cdsex3 = 3 * cdse
        self.assertEqual(12, len(cdsex3))
        self.assertEqual(3 * "Cd Cd Se Se".split(), [a.element for a in cdsex3])
        self.assertTrue(numpy.array_equal(3 * [a.xyz for a in cdse], [a.xyz for a in cdsex3]))
        self.assertFalse(set(cdse).intersection(cdsex3))
        self.assertFalse(cdse.lattice is cdsex3.lattice)
        return

    def test___imul__(self):
        """check Structure.__imul__()"""
        cdse = Structure(filename=self.cdsefile)
        lat = cdse.lattice
        els = cdse.element
        xyz = cdse.xyz
        lst = cdse.tolist()
        cdse *= 2
        self.assertEqual(8, len(cdse))
        self.assertEqual(lst, cdse[:4].tolist())
        self.assertEqual(numpy.tile(els, 2).tolist(), cdse.element.tolist())
        self.assertTrue(numpy.array_equal(numpy.tile(xyz, (2, 1)), cdse.xyz))
        self.assertEqual(8, len(set(cdse)))
        self.assertEqual(8 * [lat], [a.lattice for a in cdse])
        self.stru *= -3
        self.assertEqual(0, len(self.stru))
        return

    def test__get_lattice(self):
        """check Structure._get_lattice()"""
        lat = Lattice()
        stru = Structure()
        self.assertEqual((1, 1, 1, 90, 90, 90), stru.lattice.abcABG())
        stru2 = Structure(lattice=lat)
        self.assertTrue(lat is stru2.lattice)
        return

    def test__set_lattice(self):
        """check Structure._set_lattice()"""
        lat = Lattice()
        self.stru.lattice = lat
        self.assertEqual(2 * [lat], [a.lattice for a in self.stru])
        return

    def test_composition(self):
        """check Structure.composition property"""
        stru = self.stru
        self.assertEqual({"C": 2}, stru.composition)
        stru *= 2
        self.assertEqual({"C": 4}, stru.composition)
        stru[:] = []
        self.assertEqual({}, stru.composition)
        return

    def test_element(self):
        """check Structure.element"""
        stru = self.stru
        cdse = self.cdse
        self.assertEqual("Cd Cd Se Se".split(), cdse.element.tolist())
        self.assertEqual(cdse[:2], cdse[cdse.element == "Cd"])
        stru.element = stru.element.replace("C", "Si")
        self.assertEqual("Si", stru[0].element)
        return

    def test_xyz(self):
        """check Structure.xyz"""
        stru = self.stru
        self.assertEqual((2, 3), stru.xyz.shape)
        self.assertTrue(numpy.array_equal([1, 1, 1], stru.xyz[1]))
        stru.xyz += 0.1
        self.assertTrue(numpy.array_equal([0.1, 0.1, 0.1], stru[0].xyz))
        self.assertTrue(numpy.array_equal([1.1, 1.1, 1.1], stru[1].xyz))
        stru.xyz = 0
        stru[1].xyz[:] = 1
        self.assertTrue(numpy.array_equal([0, 0, 0], stru[0].xyz))
        self.assertTrue(numpy.array_equal([1, 1, 1], stru[1].xyz))
        # verify noop when changing empty slice
        xyz0 = numpy.copy(stru.xyz)
        stru[1:1].xyz += 1
        self.assertTrue(numpy.array_equal(xyz0, stru.xyz))
        return

    def test_x(self):
        """check Structure.x"""
        cdse = self.cdse
        self.assertEqual((4,), cdse.x.shape)
        self.assertAlmostEqual(0.6666, cdse.x[3], 5)
        stru = self.stru
        stru.x = [3, 4]
        self.assertEqual(3, stru[0].xyz[0])
        self.assertEqual(4, stru[1].xyz[0])
        return

    def test_y(self):
        """check Structure.y"""
        cdse = self.cdse
        self.assertEqual((4,), cdse.y.shape)
        self.assertAlmostEqual(0.3333, cdse.y[3], 5)
        stru = self.stru
        stru.y = [3, 4]
        self.assertEqual(3, stru[0].xyz[1])
        self.assertEqual(4, stru[1].xyz[1])
        return

    def test_z(self):
        """check Structure.z"""
        cdse = self.cdse
        self.assertEqual((4,), cdse.z.shape)
        self.assertAlmostEqual(0.87667, cdse.z[3], 5)
        stru = self.stru
        stru.z = [3, 4]
        self.assertEqual(3, stru[0].xyz[2])
        self.assertEqual(4, stru[1].xyz[2])
        return

    def test_label(self):
        """check Structure.label"""
        cdse = Structure(filename=self.cdsefile)
        self.assertEqual(4 * [""], cdse.label.tolist())
        cdse.assignUniqueLabels()
        self.assertEqual("Cd1 Cd2 Se1 Se2".split(), cdse.label.tolist())
        cdse.label = cdse.label.lower()
        self.assertEqual("cd1 cd2 se1 se2".split(), cdse.label.tolist())
        return

    def test_occupancy(self):
        """check Structure.occupancy"""
        cdse = self.cdse
        self.assertTrue(numpy.array_equal(numpy.ones(4), cdse.occupancy))
        self.stru.occupancy *= 0.5
        self.assertEqual(1.0, sum([a.occupancy for a in self.stru]))
        cdse.occupancy = 1
        self.assertTrue(all(isinstance(a.occupancy, int) for a in cdse))
        return

    def test_xyz_cartn(self):
        """check Structure.xyz_cartn"""
        pbte = copy.copy(self.pbte)
        self.assertEqual((8, 3), pbte.xyz_cartn.shape)
        self.assertTrue(numpy.allclose(6.461 / 2.0 * numpy.ones(3), pbte.xyz_cartn[0]))
        pbte.xyz_cartn += numpy.array([0.1, 0.2, 0.3]) * 6.461
        self.assertTrue(numpy.allclose([0.6, 0.7, 0.8], pbte[0].xyz))
        self.assertTrue(numpy.allclose([0.6, 0.7, 0.3], pbte[7].xyz))
        return

    def test_anisotropy(self):
        """check Structure.anisotropy"""
        self.assertEqual((2,), self.stru.anisotropy.shape)
        self.assertFalse(numpy.any(self.stru.anisotropy))
        tei = copy.copy(self.tei)
        self.assertTrue(numpy.all(tei.anisotropy))
        tei.anisotropy = False
        self.assertFalse(numpy.any(tei.anisotropy))
        self.assertAlmostEqual(0.019227, tei[0].U11, 6)
        self.assertAlmostEqual(0.019227, tei[0].U22, 6)
        self.assertAlmostEqual(0.019227, tei[0].U33, 6)
        self.assertAlmostEqual(0.0, tei[0].U12, 6)
        self.assertAlmostEqual(0.019227 * -numpy.cos(numpy.radians(128.09)), tei[0].U13, 6)
        self.assertAlmostEqual(0.0, tei[0].U23, 6)
        self.assertAlmostEqual(0.019227, tei[0].Uisoequiv, 6)
        return

    def test_U(self):
        """check Structure.U"""
        stru = self.stru
        self.assertEqual((2, 3, 3), stru.U.shape)
        self.assertFalse(numpy.any(stru.anisotropy))
        stru.U = numpy.identity(3)
        self.assertEqual(2, len(set([id(a.U) for a in stru])))
        isou = stru.lattice.isotropicunit
        self.assertTrue(numpy.array_equal(2 * [isou], stru.U))
        self.assertFalse(numpy.any(stru.anisotropy))
        stru.anisotropy = True
        stru.U = numpy.identity(3)
        self.assertTrue(numpy.array_equal(2 * [numpy.identity(3)], stru.U))
        self.assertTrue(numpy.all(stru.anisotropy))
        stru.U = 0
        self.assertTrue(numpy.all(stru.anisotropy))
        self.assertFalse(numpy.any(stru.U != 0.0))
        stru[1].U[:] = 1
        self.assertTrue(numpy.all(stru[0].U == 0.0))
        self.assertTrue(numpy.all(stru[1].U == 1.0))
        return

    def test_Uisoequiv(self):
        """check Structure.Uisoequiv"""
        tei = copy.copy(self.tei)
        self.assertEqual((16,), tei.Uisoequiv.shape)
        self.assertAlmostEqual(0.019227, tei.Uisoequiv[0], 6)
        self.assertAlmostEqual(0.019784, tei.Uisoequiv[4], 6)
        self.assertAlmostEqual(0.024813, tei.Uisoequiv[8], 6)
        self.assertAlmostEqual(0.026878, tei.Uisoequiv[12], 6)
        u11old = tei[0].U11
        tei.Uisoequiv = 0.001
        self.assertAlmostEqual(u11old * 0.001 / 0.019227, tei[0].U[0, 0])
        return

    def test_Uij(self):
        """check Structure.Uij"""
        stru = self.stru
        stru[1].anisotropy = True
        stru[1].U = [[1.1, 0.12, 0.13], [0.12, 2.2, 0.23], [0.13, 0.23, 3.3]]
        self.assertTrue(numpy.array_equal([0, 1.1], stru.U11))
        self.assertTrue(numpy.array_equal([0, 2.2], stru.U22))
        self.assertTrue(numpy.array_equal([0, 3.3], stru.U33))
        self.assertTrue(numpy.array_equal([0, 0.12], stru.U12))
        self.assertTrue(numpy.array_equal([0, 0.13], stru.U13))
        self.assertTrue(numpy.array_equal([0, 0.23], stru.U23))
        stru.U11 = stru.U22 = stru.U33 = stru.U12 = stru.U13 = stru.U23 = 0.0
        self.assertFalse(numpy.any(stru.U != 0.0))
        return

    def test_Bisoequiv(self):
        """check Structure.Bisoequiv"""
        utob = 8 * numpy.pi**2
        tei = copy.copy(self.tei)
        self.assertEqual((16,), tei.Bisoequiv.shape)
        self.assertAlmostEqual(utob * 0.019227, tei.Bisoequiv[0], 4)
        self.assertAlmostEqual(utob * 0.019784, tei.Bisoequiv[4], 4)
        self.assertAlmostEqual(utob * 0.024813, tei.Bisoequiv[8], 4)
        self.assertAlmostEqual(utob * 0.026878, tei.Bisoequiv[12], 4)
        b11old = tei[0].B11
        tei.Bisoequiv = 0.1
        self.assertAlmostEqual(b11old * 0.1 / (utob * 0.019227), tei[0].B11, 5)
        return

    def test_Bij(self):
        """check Structure.Bij"""
        stru = self.stru
        stru[1].anisotropy = True
        stru[1].U = [[1.1, 0.12, 0.13], [0.12, 2.2, 0.23], [0.13, 0.23, 3.3]]
        stru[1].U /= 8 * numpy.pi**2
        self.assertTrue(numpy.allclose([0, 1.1], stru.B11))
        self.assertTrue(numpy.allclose([0, 2.2], stru.B22))
        self.assertTrue(numpy.allclose([0, 3.3], stru.B33))
        self.assertTrue(numpy.allclose([0, 0.12], stru.B12))
        self.assertTrue(numpy.allclose([0, 0.13], stru.B13))
        self.assertTrue(numpy.allclose([0, 0.23], stru.B23))
        stru.B11 = stru.B22 = stru.B33 = stru.B12 = stru.B13 = stru.B23 = 0.0
        self.assertFalse(numpy.any(stru.U != 0.0))
        return

    def test_pickling(self):
        """Make sure Atom in Structure can be consistently pickled."""
        stru = self.stru
        a = stru[0]
        self.assertTrue(a is stru[0])
        a1, stru1 = pickle.loads(pickle.dumps((a, stru)))
        self.assertTrue(a1 is stru1[0])
        return


# End of class TestStructure

# ----------------------------------------------------------------------------

if __name__ == "__main__":
    unittest.main()

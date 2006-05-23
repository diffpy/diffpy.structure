#!/usr/bin/env python

__id__ = "$Id$"

import sys, os
import unittest
try:
    from Structure.structure import Structure, InvalidStructureFormat
    from Structure.lattice import Lattice
    from Structure.atom import Atom
except:
    print >> sys.stderr, "parent directory should be in PYTHONPATH"
    raise

##############################################################################
class TestPackageIds(unittest.TestCase):
    """test _get_package_id() in Parsers and Structure packages"""

    def setUp(self):
        self.stru = Structure()
        self.format = "cif"
        self.places = 6

    def test_Parsers_get_package_id(self):
        """check _get_package_id() in Parsers package"""
        import Structure.Parsers
        parsers_id = Structure.Parsers._get_package_id()

    def test_Structure_get_package_id(self):
        """check _get_package_id() in Structure package"""
        import Structure
        structure_id = Structure._get_package_id()

# End of TestPackageIds

##############################################################################
class TestLattice(unittest.TestCase):
    """test methods of Lattice class"""

    def setUp(self):
        self.lattice = Lattice()
        self.places = 12

    def assertListAlmostEqual(self, l1, l2, places=None):
        """wrapper for list comparison"""
        if places is None: places = self.places
        self.assertEqual(len(l1), len(l2))
        for i in range(len(l1)):
            self.assertAlmostEqual(l1[i], l2[i], places)

    def test_setLatPar(self):
        """check calculation of standard unit cell vectors"""
        from numarray import dot
        from math import radians, sqrt, cos, sin
        norm = lambda x : sqrt(sum([xi**2 for xi in x]))
        cosd = lambda x : cos(radians(x))
        sind = lambda x : sin(radians(x))
        self.lattice.setLatPar(1.0, 2.0, 3.0, 80, 100, 120)
        base = self.lattice.base
        self.assertAlmostEqual(1.0, norm(base[0]), self.places)
        self.assertAlmostEqual(2.0, norm(base[1]), self.places)
        self.assertAlmostEqual(3.0, norm(base[2]), self.places)
        self.assertAlmostEqual(cosd(80.0),
                dot(base[1],base[2])/(2*3), self.places)
        self.assertAlmostEqual(cosd(100.0),
                dot(base[0],base[2])/(1*3), self.places)
        self.assertAlmostEqual(cosd(120.0),
                dot(base[0],base[1])/(1*2), self.places)

    def test_setLatBase(self):
        """check calculation of unit cell rotation"""
        import numarray as num
        import numarray.linear_algebra as numalg
        base = num.array([[ 1.0,  1.0,  0.0],
                          [ 0.0,  1.0,  1.0],
                          [ 1.0,  0.0,  1.0]])
        self.lattice.setLatBase(base)
        self.assertAlmostEqual(self.lattice.a, num.sqrt(2.0), self.places)
        self.assertAlmostEqual(self.lattice.b, num.sqrt(2.0), self.places)
        self.assertAlmostEqual(self.lattice.c, num.sqrt(2.0), self.places)
        self.assertAlmostEqual(self.lattice.alpha, 60.0, self.places)
        self.assertAlmostEqual(self.lattice.beta,  60.0, self.places)
        self.assertAlmostEqual(self.lattice.gamma, 60.0, self.places)
        detR0 = numalg.determinant(self.lattice.baserot)
        self.assertAlmostEqual(detR0, 1.0, self.places)
        # try if rotation matrix works
        self.assertEqual(num.all(base == self.lattice.base), True)
        self.lattice.setLatPar(alpha=44, beta=66, gamma=88)
        self.assertNotEqual(num.all(base == self.lattice.base), True)
        self.lattice.setLatPar(alpha=60, beta=60, gamma=60)
        self.assertListAlmostEqual(base[0], self.lattice.base[0])
        self.assertListAlmostEqual(base[1], self.lattice.base[1])
        self.assertListAlmostEqual(base[2], self.lattice.base[2])

    def test_repr(self):
        """check string representation of this lattice"""
        r = repr(self.lattice)
        self.assertEqual(r, "Lattice()")
        self.lattice.setLatPar(1, 2, 3, 10, 20, 30)
        r = repr(self.lattice)
        r0 = "Lattice(a=1, b=2, c=3, alpha=10, beta=20, gamma=30)"
        self.assertEqual(r, r0)
        base = [[ 1.0,  1.0,  0.0],
                [ 0.0,  2.0,  2.0],
                [ 3.0,  0.0,  3.0]]
        self.lattice.setLatBase(base)
        r = repr(self.lattice)
        self.assertEqual(r, "Lattice(base=%r)" % self.lattice.base)

# End of TestLattice

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

    def test_cartesian(self):
        """check conversion to cartesian coordinates"""
        from math import sqrt
        stru = self.stru
        s_rc0 = stru.cartesian(stru[0])
        f_rc0 = 3*[0.0]
        s_rc1 = stru.cartesian(stru[1])
        f_rc1 = [sqrt(0.75), 0.5, 1.]
        self.assertListAlmostEqual(s_rc0, f_rc0)
        self.assertListAlmostEqual(s_rc1, f_rc1)

    def test_placeInLattice(self):
        """check conversion of coordinates when placed in new lattice"""
        stru = self.stru
        new_lattice = Lattice(.5, .5, .5, 90, 90, 60)
        stru.placeInLattice(new_lattice)
        a0 = stru[0]
        self.assertListAlmostEqual(a0.xyz, 3*[0.0])
        a1 = stru[1]
        self.assertListAlmostEqual(a1.xyz, [2.0, 0.0, 2.0])
# End of TestStructure

##############################################################################
class TestP_xyz(unittest.TestCase):
    """test Parser for xyz file format"""

    def setUp(self):
        self.stru = Structure()
        self.format = 'xyz'
        import os, tempfile
        handle, self.tmpname = tempfile.mkstemp()
        os.close(handle)

    def tearDown(self):
        import os
        os.remove(self.tmpname)

    def test_read_xyz(self):
        """check reading of normal xyz file"""
        stru = self.stru
        stru.read('testdata/bucky.xyz', self.format)
        s_els = [a.element for a in stru]
        self.assertEqual(stru.title, 'bucky-ball')
        self.assertEqual(s_els, 60*['C'])

    def test_read_xyz_bad(self):
        """check exceptions when reading invalid xyz file"""
        stru = self.stru
        self.assertRaises(InvalidStructureFormat, stru.read,
                'testdata/bucky-bad1.xyz', self.format )
        self.assertRaises(InvalidStructureFormat, stru.read,
                'testdata/bucky-bad2.xyz', self.format )
        self.assertRaises(InvalidStructureFormat, stru.read,
                'testdata/bucky-plain.xyz', self.format )
        self.assertRaises(InvalidStructureFormat, stru.read,
                'testdata/hexagon-raw.xy', self.format )

    def test_writeStr_xyz(self):
        """check string representation of normal xyz file"""
        import re
        stru = self.stru
        stru.title = "test of writeStr"
        stru.lattice = Lattice(1.0, 2.0, 3.0, 90.0, 90.0, 90.0)
        stru[:] = [
            Atom('H', [1., 1., 1.]),
            Atom('Cl', [3., 2., 1.])
        ]
        s1 = stru.writeStr(self.format)
        s1 = re.sub('[ \t]+', ' ', s1)
        s0 = "2\n%s\nH 1 2 3\nCl 3 4 3\n" % stru.title
        self.assertEqual(s1, s0)

    def test_write_xyz(self):
        """check writing of normal xyz file"""
        import re
        stru = self.stru
        stru.title = "test of writeStr"
        stru.lattice = Lattice(1.0, 2.0, 3.0, 90.0, 90.0, 90.0)
        stru[:] = [
            Atom('H', [1., 1., 1.]),
            Atom('Cl', [3., 2., 1.])
        ]
        stru.write(self.tmpname, self.format)
        f_s = open(self.tmpname).read()
        f_s = re.sub('[ \t]+', ' ', f_s)
        s_s = "2\n%s\nH 1 2 3\nCl 3 4 3\n" % stru.title
        self.assertEqual(f_s, s_s)
# End of TestP_xyz

##############################################################################
class TestP_rawxyz(unittest.TestCase):
    """test Parser for rawxyz file format"""

    def setUp(self):
        self.stru = Structure()
        self.format = "rawxyz"

    def test_read_plainxyz(self):
        """check reading of a plain xyz file"""
        stru = self.stru
        stru.read('testdata/bucky-plain.xyz', self.format)
        s_els = [a.element for a in stru]
        self.assertEqual(stru.title, '')
        self.assertEqual(s_els, 60*['C'])

    def test_read_plainxyz_bad(self):
        """check exceptions when reading invalid plain xyz file"""
        stru = self.stru
        self.assertRaises(InvalidStructureFormat, stru.read,
                'testdata/bucky-plain-bad.xyz', self.format)

    def test_read_rawxyz(self):
        """check reading of raw xyz file"""
        stru = self.stru
        stru.read('testdata/bucky-raw.xyz', self.format)
        s_els = [a.element for a in stru]
        self.assertEqual(stru.title, "")
        self.assertEqual(s_els, 60*[''])
        stru.read('testdata/hexagon-raw.xyz', self.format)
        zs = [a.xyz[-1] for a in stru]
        self.assertEqual(zs, 6*[0.0])

    def test_read_rawxyz_bad(self):
        """check exceptions when reading unsupported xy file"""
        stru = self.stru
        self.assertRaises(InvalidStructureFormat, stru.read,
                'testdata/hexagon-raw-bad.xyz', self.format)
        self.assertRaises(InvalidStructureFormat, stru.read,
                'testdata/hexagon-raw.xy', self.format)

    def test_writeStr_rawxyz(self):
        """check writing of normal xyz file"""
        import re
        stru = self.stru
        stru.title = "test of writeStr"
        stru.lattice = Lattice(1.0, 2.0, 3.0, 90.0, 90.0, 90.0)
        # plain version
        stru[:] = [ Atom('H', [1., 1., 1.]) ]
        s1 = stru.writeStr(self.format)
        s1 = re.sub('[ \t]+', ' ', s1)
        s0 = "H 1 2 3\n"
        # brutal raw version
        stru[0].element = ""
        s1 = stru.writeStr(self.format)
        s0 = "1 2 3\n"
        self.assertEqual(s1, s0)
# End of TestP_rawxyz

##############################################################################
class TestP_pdffit(unittest.TestCase):
    """test Parser for PDFFit file format"""

    def setUp(self):
        self.stru = Structure()
        self.format = "pdffit"
        self.places = 8

    def assertListAlmostEqual(self, l1, l2, places=None):
        """wrapper for list comparison"""
        if places is None: places = self.places
        self.assertEqual(len(l1), len(l2))
        for i in range(len(l1)):
            self.assertAlmostEqual(l1[i], l2[i], places)

    def test_read_pdffit_ZnSb(self):
        """check reading of ZnSb pdffit structure file"""
        stru = self.stru
        stru.read('testdata/ZnSb_RT_Q28X_VM_20_fxiso.rstr', self.format)
        f_title = "Cell structure file of Zn4Sb3 #167 interstitial"
        self.assertEqual(stru.title, f_title)
        self.assertAlmostEqual(stru.pdffit['scale'], 0.826145)
        self.assertAlmostEqual(stru.pdffit['delta'], 4.687951)
        self.assertAlmostEqual(stru.pdffit['gamma'], 0.01)
        self.assertAlmostEqual(stru.pdffit['srat'], 1.02)
        self.assertAlmostEqual(stru.pdffit['rcut'], 0.03)
        self.assertEqual(stru.pdffit['spcgr'], 'R-3c')
        s_lat = [ stru.lattice.a, stru.lattice.b, stru.lattice.c,
            stru.lattice.alpha, stru.lattice.beta, stru.lattice.gamma ]
        f_lat = [12.309436, 12.309436, 12.392839, 90.0, 90.0, 120.0]
        self.assertListAlmostEqual(s_lat, f_lat)
        s_dcell = stru.pdffit['dcell']
        f_dcell = [0.000008, 0.000008, 0.000013, 0.0, 0.0, 0.0]
        self.assertListAlmostEqual(s_dcell, f_dcell)
        self.assertEqual(stru.pdffit['ncell'], [1,1,1,66])
        s_els = [a.element for a in stru]
        self.assertEqual(s_els, 36*['Zn']+30*['Sb'])
        a0 = stru[0]
        s_xyz = a0.xyz
        f_xyz = [0.09094387, 0.24639539, 0.40080261];
        s_o = a0.occupancy
        f_o = 0.9
        s_sigxyz = a0.sigxyz
        f_sigxyz = [ 0.00000079, 0.00000076, 0.00000064];
        s_sigo = a0.sigo
        f_sigo = 0.0
        s_U = [ a0.U[i][i] for i in range(3) ]
        f_U = 3*[0.01]
        self.assertListAlmostEqual(s_xyz, f_xyz)
        self.assertListAlmostEqual(s_sigxyz, f_sigxyz)
        self.assertListAlmostEqual(s_U, f_U)
        self.assertAlmostEqual(s_o, f_o)
        self.assertAlmostEqual(s_sigo, f_sigo)

    def test_read_pdffit_Ni(self):
        """check reading of Ni pdffit structure file"""
        stru = self.stru
        stru.read('testdata/Ni.stru', self.format)
        f_title = "structure Ni  FCC"
        self.assertEqual(stru.title, f_title)
        self.assertEqual(stru.pdffit['spcgr'], 'Fm-3m')
        s_lat = [ stru.lattice.a, stru.lattice.b, stru.lattice.c,
            stru.lattice.alpha, stru.lattice.beta, stru.lattice.gamma ]
        f_lat = [3.52, 3.52, 3.52, 90.0, 90.0, 90.0]
        for i in range(len(s_lat)):
            self.assertAlmostEqual(s_lat[i], f_lat[i])
        self.assertEqual(stru.pdffit['ncell'], [1,1,1,4])
        s_els = [a.element for a in stru]
        self.assertEqual(s_els, 4*['Ni'])
        a0 = stru[0]
        s_xyz = a0.xyz
        f_xyz = [0.0, 0.0, 0.0];
        s_o = a0.occupancy
        f_o = 1.0
        s_U = [ a0.U[i][i] for i in range(3) ]
        f_U = 3*[0.00126651]
        for i in range(3):
            self.assertAlmostEqual(s_xyz[i], f_xyz[i])
            self.assertAlmostEqual(s_U[i], f_U[i])
        self.assertAlmostEqual(s_o, f_o)

    def test_read_pdffit_Ni_prim123(self):
        """check reading of Ni_prim supercell 1x2x3"""
        stru = self.stru
        stru.read('testdata/Ni_prim123.stru', self.format)
        s_lat = [ stru.lattice.a, stru.lattice.b, stru.lattice.c,
            stru.lattice.alpha, stru.lattice.beta, stru.lattice.gamma ]
        f_lat = [2.489016,  2*2.489016,  3*2.489016, 60.0, 60.0, 60.0]
        for i in range(len(s_lat)):
            self.assertAlmostEqual(s_lat[i], f_lat[i])
        s_els = [a.element for a in stru]
        self.assertEqual(s_els, 6*['Ni'])
        a5 = stru[5]
        s_xyz = a5.xyz
        f_xyz = [0.0, 1.0/2.0, 2.0/3.0];
        for i in range(3):
            self.assertAlmostEqual(s_xyz[i], f_xyz[i])
        s_o = a5.occupancy
        f_o = 1.0
        self.assertAlmostEqual(s_o, f_o)
        s_U = [ a5.U[ij[0],ij[1]]
                for ij in [(0,0), (1,1), (2,2), (0,1), (0,2), (1,2)] ]
        f_U = 3*[0.00126651] + 3*[-0.00042217]
        for i in range(len(s_U)):
            self.assertAlmostEqual(s_U[i], f_U[i])

    def test_read_pdffit_bad(self):
        """check exceptions when reading invalid pdffit file"""
        stru = self.stru
        self.assertRaises(InvalidStructureFormat, stru.read,
                'testdata/Ni-bad.stru', self.format)
        self.assertRaises(InvalidStructureFormat, stru.read,
                'testdata/bucky.xyz', self.format)

    def test_writeStr_pdffit(self):
        """check writing of normal xyz file"""
        import re
        stru = self.stru
        stru.read('testdata/Ni.stru', self.format)
        f_s = open('testdata/Ni.stru').read()
        f_s = re.sub('[ \t]+', ' ', f_s)
        f_s = re.sub('[ \t]+\n', '\n', f_s)
        s_s = stru.writeStr(self.format)
        s_s = re.sub('[ \t]+', ' ', s_s)
        self.assertEqual(f_s, s_s)
# End of TestP_pdffit

##############################################################################
class TestP_pdb(unittest.TestCase):
    """test Parser for PDB file format"""

    def setUp(self):
        self.stru = Structure()
        self.format = "pdb"
        self.places = 3

    def assertListAlmostEqual(self, l1, l2, places=None):
        """wrapper for list comparison"""
        if places is None: places = self.places
        self.assertEqual(len(l1), len(l2))
        for i in range(len(l1)):
            self.assertAlmostEqual(l1[i], l2[i], places)

    def test_read_pdb_arginine(self):
        """check reading of arginine PDB file"""
        stru = self.stru
        stru.read('testdata/arginine.pdb', self.format)
        f_els = [ "N", "C", "C", "O", "C", "C", "C", "N", "C", "N", "N", "H",
            "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H",
            "O", "H" ]
        s_els = [a.element for a in stru]
        self.assertEqual(s_els, f_els)
        s_lat = [ stru.lattice.a, stru.lattice.b, stru.lattice.c,
            stru.lattice.alpha, stru.lattice.beta, stru.lattice.gamma ]
        f_lat = [1.0, 1.0, 1.0, 90.0, 90.0, 90.0]
        self.assertEqual(s_lat, f_lat)
        a0 = stru[0]
        self.assertListAlmostEqual(a0.xyz, [0.735, 2.219, 1.389])

    def test_rwStr_pdb_CdSe(self):
        """check conversion to PDB file format"""
        stru = self.stru
        stru.read('testdata/CdSe_bulk.stru', 'pdffit')
        s = stru.writeStr(self.format)
        # all lines should be 80 characters long
        linelens = [ len(l) for l in s.split('\n') if l != "" ]
        self.assertEqual(linelens, len(linelens)*[80])
        # now clean and re-read structure
        stru = Structure()
        stru.readStr(s, self.format)
        s_els = [a.element for a in stru]
        f_els = ['Cd', 'Cd', 'Se', 'Se']
        self.assertEqual(s_els, f_els)
        s_lat = [ stru.lattice.a, stru.lattice.b, stru.lattice.c,
            stru.lattice.alpha, stru.lattice.beta, stru.lattice.gamma ]
        f_lat = [ 4.235204,  4.235204,  6.906027, 90.0, 90.0, 120.0 ]
        self.assertListAlmostEqual(s_lat, f_lat)
        a0 = stru[0]
        s_Uii = [ a0.U[i,i] for i in range(3) ]
        f_Uii = [ 0.01303035, 0.01303035, 0.01401959 ]
        self.assertListAlmostEqual(s_Uii, f_Uii)
        s_sigUii = [ a0.sigU[i,i] for i in range(3) ]
        f_sigUii = [ 0.00011127, 0.00011127, 0.00019575 ]
        self.assertListAlmostEqual(s_sigUii, f_sigUii)
        s_title = stru.title
        f_title = "Cell structure file of CdSe #186"
        self.assertEqual(s_title, f_title)

# End of TestP_pdb

##############################################################################
class TestP_xcfg(unittest.TestCase):
    """test Parser for XCFG file format"""

    def setUp(self):
        self.stru = Structure()
        self.format = "xcfg"
        self.places = 6

    def assertListAlmostEqual(self, l1, l2, places=None):
        """wrapper for list comparison"""
        if places is None: places = self.places
        self.assertEqual(len(l1), len(l2))
        for i in range(len(l1)):
            self.assertAlmostEqual(l1[i], l2[i], places)

    def test_read_xcfg(self):
        """check reading of BubbleRaft XCFG file"""
        stru = self.stru
        stru.read('testdata/BubbleRaftShort.xcfg', self.format)
        f_els = 500* [ "Ar" ]
        s_els = [a.element for a in stru]
        self.assertEqual(s_els, f_els)
        self.assertAlmostEqual(stru.dist(stru[82],stru[357]), 47.5627, 3)
        s_lat = [ stru.lattice.a, stru.lattice.b, stru.lattice.c,
            stru.lattice.alpha, stru.lattice.beta, stru.lattice.gamma ]
        f_lat = [127.5, 119.5, 3.0, 90.0, 90.0, 90.0]
        self.assertListAlmostEqual(s_lat, f_lat)

    def test_rwStr_xcfg_CdSe(self):
        """check conversion to XCFG file format"""
        stru = self.stru
        stru.read('testdata/CdSe_bulk.stru', 'pdffit')
        s = stru.writeStr(self.format)
        stru = Structure()
        stru.readStr(s, self.format)
        s_els = [a.element for a in stru]
        f_els = ['Cd', 'Cd', 'Se', 'Se']
        self.assertEqual(s_els, f_els)
        s_lat = [ stru.lattice.a, stru.lattice.b, stru.lattice.c,
            stru.lattice.alpha, stru.lattice.beta, stru.lattice.gamma ]
        f_lat = [ 4.235204,  4.235204,  6.906027, 90.0, 90.0, 120.0 ]
        self.assertListAlmostEqual(s_lat, f_lat)
        a0 = stru[0]
        s_Uii = [ a0.U[i,i] for i in range(3) ]
        f_Uii = [ 0.01303035, 0.01303035, 0.01401959 ]
        self.assertListAlmostEqual(s_Uii, f_Uii)

# End of TestP_xcfg

##############################################################################
class TestP_cif(unittest.TestCase):
    """test Parser for CIF file format"""

    def setUp(self):
        self.stru = Structure()
        self.format = "cif"
        self.places = 6

    def test_writeStr_cif(self):
        """check conversion to CIF string"""
        stru = self.stru
        stru.read('testdata/CdSe_bulk.stru', 'pdffit')
        s_s = stru.writeStr(self.format)

# End of TestP_cif

##############################################################################
def ParsersSuite():
    suite = unittest.TestSuite()
    suite.addTest( unittest.makeSuite(TestLattice) )
    suite.addTest( unittest.makeSuite(TestStructure) )
    suite.addTest( unittest.makeSuite(TestP_xyz) )
    suite.addTest( unittest.makeSuite(TestP_rawxyz) )
    suite.addTest( unittest.makeSuite(TestP_pdffit) )
    suite.addTest( unittest.makeSuite(TestP_pdb) )
    suite.addTest( unittest.makeSuite(TestP_xcfg) )
    suite.addTest( unittest.makeSuite(TestP_cif) )
    return suite

if __name__ == '__main__':
    # change to directory of this script
    os.chdir(os.path.abspath(os.path.dirname(sys.argv[0])))
#   seltests = unittest.TestSuite()
#   seltests.addTest( ParsersSuite() )
#   seltests.addTest( unittest.makeSuite(TestLattice) )
    # this runs all tests:
    unittest.main()
    # and this only some:
#   unittest.TextTestRunner().run(seltests)
#   unittest.TextTestRunner(verbosity=100).run(seltests)

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

from diffpy.Structure import Structure, InvalidStructureFormat
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

if __name__ == '__main__':
    unittest.main()

# End of file

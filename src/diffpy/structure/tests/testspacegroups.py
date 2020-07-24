#!/usr/bin/env python
##############################################################################
#
# diffpy.structure  Complex Modeling Initiative
#                   (c) 2019 Brookhaven Science Associates,
#                   Brookhaven National Laboratory.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################

"""Unit tests for diffpy.structure.spacegroups
"""


import unittest

from diffpy.structure.spacegroups import SpaceGroupList, _hashSymOpList
from diffpy.structure.spacegroups import GetSpaceGroup, FindSpaceGroup

# ----------------------------------------------------------------------------

class TestRoutines(unittest.TestCase):

    def setUp(self):
        return


    def test_old_alt_name(self):
        "check GetSpaceGroup lookup from deprecated alt_name values"
        # alt_name values which do not map to short_name or pdb_name
        altnames_sgnos = (
            ("P M 3", 200), ("P N 3", 201), ("F M 3", 202),
            ("F D 3", 203), ("I M 3", 204), ("P A 3", 205),
            ("I A 3", 206), ("P M 3 M", 221), ("P N 3 N", 222),
            ("P M 3 N", 223), ("P N 3 M", 224), ("F M 3 M", 225),
            ("F M 3 C", 226), ("F D 3 M", 227), ("F D 3 C", 228),
            ("I M 3 M", 229), ("I A 3 D", 230)
        )
        for name, sgno in altnames_sgnos:
            self.assertIs(GetSpaceGroup(sgno), GetSpaceGroup(name))
        return


    def test_GetSpaceGroup(self):
        "check GetSpaceGroup function"
        from diffpy.structure.spacegroups import sg125
        self.assertRaises(ValueError, GetSpaceGroup, 0)
        self.assertRaises(ValueError, GetSpaceGroup, 300)
        self.assertRaises(ValueError, GetSpaceGroup, "300")
        self.assertIs(sg125, GetSpaceGroup(125))
        self.assertIs(sg125, GetSpaceGroup("125"))
        self.assertIs(sg125, GetSpaceGroup("P4/nbm"))
        self.assertIs(sg125, GetSpaceGroup("P 4/n 2/b 2/m"))
        # old alt_name
        self.assertIs(sg125, GetSpaceGroup("P 4/N B M"))
        # upper case pdb_name
        self.assertIs(sg125, GetSpaceGroup("P 4/N 2/B 2/M"))
        return


    def test_FindSpaceGroup(self):
        "check FindSpaceGroup function"
        sg123 = GetSpaceGroup(123)
        ops123 = list(sg123.iter_symops())
        self.assertRaises(ValueError, FindSpaceGroup, [])
        self.assertRaises(ValueError, FindSpaceGroup, 2 * ops123)
        self.assertIs(sg123, FindSpaceGroup(ops123))
        sg123r = FindSpaceGroup(ops123[::-1])
        self.assertIsNot(sg123, sg123r)
        self.assertIsNot(sg123.symop_list, sg123r.symop_list)
        self.assertEqual(ops123[::-1], sg123r.symop_list)
        self.assertEqual(_hashSymOpList(sg123.symop_list),
                         _hashSymOpList(sg123r.symop_list))
        self.assertIs(sg123, FindSpaceGroup(ops123[::-1], shuffle=True))
        return


    def test__hashSymOpList(self):
        "verify _hashSymOpList is unique for each spacegroup"
        hset = set(_hashSymOpList(sg.symop_list) for sg in SpaceGroupList)
        self.assertEqual(len(SpaceGroupList), len(hset))
        return

    def test_spacegroup_representation(self):
        """Verify SpaceGroup.__repr__()."""
        numbers = [1, 3, 16, 75, 143, 168, 229]
        short_names = ["P1", "P2", "P222", "P4", "P3", "P6", "Im-3m"]
        systems = [
            "Triclinic",
            "Monoclinic",
            "Orthorhombic",
            "Tetragonal",
            "Trigonal",
            "Hexagonal",
            "Cubic",
        ]
        n_symmetry_matrices = [1, 2, 4, 4, 3, 6, 96]
        n_point_symmetry_matrices = [1, 2, 4, 4, 3, 6, 48]
        sg_dict = dict(zip(
            numbers,
            zip(short_names, systems, n_symmetry_matrices, n_point_symmetry_matrices)
        ))
        for key, value in sg_dict.items():
            sg = GetSpaceGroup(key)
            expected_str = (
                               "SpaceGroup #%i (%s, %s). Symmetry matrices: %i, "
                               "point sym. matr.: %i"
                           ) % ((key,) + value)
            self.assertEqual(sg.__repr__(), expected_str)
        return

# End of class TestRoutines

# ----------------------------------------------------------------------------

if __name__ == '__main__':
    unittest.main()

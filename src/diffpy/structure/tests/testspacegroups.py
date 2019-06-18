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

# End of class TestRoutines

# ----------------------------------------------------------------------------

if __name__ == '__main__':
    unittest.main()

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

# ----------------------------------------------------------------------------

class TestRoutines(unittest.TestCase):

    def setUp(self):
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

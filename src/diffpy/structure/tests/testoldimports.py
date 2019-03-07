#!/usr/bin/env python
##############################################################################
#
# diffpy.structure  Complex Modeling Initiative
#                   (c) 2017 Brookhaven Science Associates,
#                   Brookhaven National Laboratory.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################

"""
Unit tests for imports of old camel-case names.
"""


import sys
import warnings
import importlib
import unittest

import diffpy

# ----------------------------------------------------------------------------

class TestOldImports(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        "Uncache any already-imported old modules."
        for modname in tuple(sys.modules):
            if modname.startswith('diffpy.Structure'):
                del sys.modules[modname]    # pragma: no cover
        return


    def test_00TopImport(self):
        """check import of diffpy.Structure
        """
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', category=DeprecationWarning)
            import diffpy.Structure as m0
        self.assertIs(diffpy.structure, m0)
        # second import should raise no warning
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            import diffpy.Structure as m1
        self.assertIs(diffpy.structure, m1)
        return


    def test_O1SubmoduleImport(self):
        """check import of diffpy.Structure submodules.
        """
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always', category=DeprecationWarning)
            import diffpy.Structure.SymmetryUtilities as symutil
            self.assertIs(DeprecationWarning, w[0].category)
        self.assertIs(diffpy.structure.symmetryutilities, symutil)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always', category=DeprecationWarning)
            import diffpy.Structure.Parsers.P_cif as pcif
            self.assertIs(DeprecationWarning, w[0].category)
        self.assertIs(diffpy.structure.parsers.p_cif, pcif)
        self.assertRaises(ImportError, importlib.import_module,
                          'diffpy.Structure.SSpaceGroups')
        return

# End of class TestOldImports

# ----------------------------------------------------------------------------

if __name__ == '__main__':
    unittest.main()

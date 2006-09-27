"""Unit tests for version-info functions in Structure module.
"""

__id__ = "$Id$"

import unittest

##############################################################################
class TestStructureModule(unittest.TestCase):
    """test _get_package_id() in Parsers and Structure packages"""

    def setUp(self):
        return

    def test_Structure__get_package_id(self):
        """check _get_package_id() in Structure package"""
        import Structure
        structure_id = Structure._get_package_id()

    def test_Parsers__get_package_id(self):
        """check _get_package_id() in Parsers package"""
        import Structure.Parsers
        parsers_id = Structure.Parsers._get_package_id()

# End of TestStructureModule

if __name__ == '__main__':
    unittest.main()

# End of file

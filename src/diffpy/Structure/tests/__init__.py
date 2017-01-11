#!/usr/bin/env python
##############################################################################
#
# diffpy.Structure  by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2010 Trustees of the Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################

"""Unit tests for diffpy.Structure.
"""

def testsuite():
    '''Build a unit tests suite for the diffpy.Structure package.

    Return a unittest.TestSuite object.
    '''
    import unittest
    modulenames = '''
        diffpy.Structure.tests.TestAtom
        diffpy.Structure.tests.TestLattice
        diffpy.Structure.tests.TestLoadStructure
        diffpy.Structure.tests.TestP_cif
        diffpy.Structure.tests.TestP_discus
        diffpy.Structure.tests.TestP_pdffit
        diffpy.Structure.tests.TestParsers
        diffpy.Structure.tests.TestStructure
        diffpy.Structure.tests.TestSuperCell
        diffpy.Structure.tests.TestSymmetryUtilities
    '''.split()
    suite = unittest.TestSuite()
    loader = unittest.defaultTestLoader
    mobj = None
    for mname in modulenames:
        exec ('import %s as mobj' % mname)
        suite.addTests(loader.loadTestsFromModule(mobj))
    return suite


def test():
    '''Execute all unit tests for the diffpy.Structure package.
    Return a unittest TestResult object.
    '''
    import unittest
    suite = testsuite()
    runner = unittest.TextTestRunner()
    result = runner.run(suite)
    return result


# End of file

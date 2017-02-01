#!/usr/bin/env python3
##############################################################################
#
# diffpy.structure  by DANSE Diffraction group
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

"""Unit tests for diffpy.structure.
"""

def testsuite():
    '''Build a unit tests suite for the diffpy.structure package.

    Return a unittest.TestSuite object.
    '''
    import unittest
    modulenames = '''
        diffpy.structure.tests.testatom
        diffpy.structure.tests.testlattice
        diffpy.structure.tests.testloadstructure
        diffpy.structure.tests.testp_cif
        diffpy.structure.tests.testp_discus
        diffpy.structure.tests.testp_pdffit
        diffpy.structure.tests.testparsers
        diffpy.structure.tests.teststructure
        diffpy.structure.tests.testsupercell
        diffpy.structure.tests.testsymmetryutilities
    '''.split()
    suite = unittest.TestSuite()
    loader = unittest.defaultTestLoader
    for mname in modulenames:
        ns = {}
        exec('import {} as mobj'.format(mname), ns)
        suite.addTests(loader.loadTestsFromModule(ns['mobj']))
    return suite


def test():
    '''Execute all unit tests for the diffpy.structure package.
    Return a unittest TestResult object.
    '''
    import unittest
    suite = testsuite()
    runner = unittest.TextTestRunner()
    result = runner.run(suite)
    return result

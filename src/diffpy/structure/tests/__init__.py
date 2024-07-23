#!/usr/bin/env python
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

import unittest


def testsuite(pattern=""):
    """Create a unit tests suite for diffpy.structure package.

    Parameters
    ----------
    pattern : str, optional
        Regular expression pattern for selecting test cases.
        Select all tests when empty.  Ignore the pattern when
        any of unit test modules fails to import.

    Returns
    -------
    suite : `unittest.TestSuite`
        The TestSuite object containing the matching tests.
    """
    import re
    from itertools import chain
    from os.path import dirname

    # Use importlib.resources for Python 3.9 and newer or importlib_resources for older versions
    try:
        from importlib import resources as importlib_resources
    except ImportError:
        import importlib_resources

    loader = unittest.defaultTestLoader
    # Use importlib_resources.files() to get a reference to the package directory
    thisdir = str(importlib_resources.files(__name__))
    depth = __name__.count(".") + 1
    topdir = thisdir
    for i in range(depth):
        topdir = dirname(topdir)
    suite_all = loader.discover(thisdir, top_level_dir=topdir)
    # always filter the suite by pattern to test-cover the selection code.
    suite = unittest.TestSuite()
    rx = re.compile(pattern)
    tsuites = list(chain.from_iterable(suite_all))
    tsok = all(isinstance(ts, unittest.TestSuite) for ts in tsuites)
    if not tsok:  # pragma: no cover
        return suite_all
    tcases = chain.from_iterable(tsuites)
    for tc in tcases:
        tcwords = tc.id().split(".")
        shortname = ".".join(tcwords[-3:])
        if rx.search(shortname):
            suite.addTest(tc)
    # verify all tests are found for an empty pattern.
    assert pattern or suite_all.countTestCases() == suite.countTestCases()
    return suite


def test():
    """Execute all unit tests for the diffpy.structure package.

    Returns
    -------
    result : `unittest.TestResult`
    """
    suite = testsuite()
    runner = unittest.TextTestRunner()
    result = runner.run(suite)
    return result

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

"""Helper routines for running other unit tests.
"""

# helper functions

def datafile(filename):
    from importlib.resources import as_file, files

    ref = files(__package__) / ("testdata/" + filename)
    with as_file(ref) as rv:
        return rv
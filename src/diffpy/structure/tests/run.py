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

"""Convenience module for executing all unit tests with

python -m diffpy.structure.tests.run
"""


if __name__ == '__main__':
    import sys
    from diffpy.structure.tests import test
    # produce zero exit code for a successful test
    sys.exit(not test().wasSuccessful())

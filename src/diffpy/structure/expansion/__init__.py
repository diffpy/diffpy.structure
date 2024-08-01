#!/usr/bin/env python
##############################################################################
#
# diffpy.structure  by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2008 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Chris Farrow, Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################

"""Methods and classes for manipulating `Structure` instances.

Package content:
    * supercell -- create a supercell from an existing `Structure`.
"""

# Import below whatever should be available at package namespace.

from diffpy.structure.expansion.supercell_mod import supercell

# silence pyflakes checker
assert supercell

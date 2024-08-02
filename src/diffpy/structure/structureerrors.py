#!/usr/bin/env python
##############################################################################
#
# diffpy.structure  by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2006 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################

"""Exceptions used in Structure package.
"""


class StructureFormatError(Exception):
    """Exception for failed IO from Structure file."""

    pass


class LatticeError(Exception):
    """Exception for impossible lattice parameters."""

    pass


class SymmetryError(Exception):
    """Exception raised for invalid symmetry operations."""

    pass

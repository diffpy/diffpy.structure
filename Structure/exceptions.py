########################################################################
#
# <PackageName>     by DANSE Diffraction group
#                   Simon J.L. Billinge
#                   Michigan State University
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See COPYRIGHT.txt for copying and usage conditions.
# See LICENSE.txt for license information.
#
########################################################################

"""Exceptions used in Structure package:
    InvalidLattice
"""

__id__ = "$Id$"

class InvalidStructureFormat(Exception):
    """Exception for failed IO from Structure file
    """

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value

# End of InvalidStructureFormat


class InvalidLattice(Exception):
    """Exception for impossible lattice parameters.
    """

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value

# End of InvalidLattice


class SymmetryError(Exception):
    """Exception raised for invalid symmetry operations.
    """

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value

# End of SymmetryError

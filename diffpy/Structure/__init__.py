#!/usr/bin/env python
##############################################################################
#
# Structure         by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2006 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################

"""classes related to structure of materials
Classes:
    Atom
    Lattice
    Structure
    PDFFitStructure
Exceptions:
    StructureFormatError
    LatticeError
    SymmetryError
"""

__id__ = "$Id$"

##############################################################################
# interface definitions
##############################################################################

from diffpy.Structure.StructureErrors import *
from diffpy.Structure.atom import Atom
from diffpy.Structure.lattice import Lattice
from diffpy.Structure.structure import Structure
from diffpy.Structure.pdffitstructure import PDFFitStructure

# obtain version information
from diffpy.Structure.version import __version__

# unit tests
def test():
    '''Execute all unit tests for the diffpy.Structure package.
    Return a unittest TestResult object.
    '''
    from diffpy.Structure.tests import test
    return test()

# End of file

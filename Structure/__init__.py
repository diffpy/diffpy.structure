########################################################################
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
########################################################################

"""classes related to structure of materials
Classes:
    Atom
    Lattice
    Structure
    PDFFitStructure
Exceptions:
    InvalidStructureFormat
    InvalidLattice
    SymmetryError
"""

__id__ = "$Id$"

##############################################################################
# interface definitions
##############################################################################

from Structure import Structure
from Lattice import Lattice
from Atom import Atom
from PDFFitStructure import PDFFitStructure
from StructureErrors import InvalidStructureFormat
from StructureErrors import InvalidLattice
from StructureErrors import SymmetryError
from version import __version__

# End of file

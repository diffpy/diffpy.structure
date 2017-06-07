#!/usr/bin/env python
##############################################################################
#
# diffpy.Structure  by DANSE Diffraction group
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

##############################################################################
# interface definitions
##############################################################################

from diffpy.Structure.StructureErrors import (
    StructureFormatError, LatticeError, SymmetryError)
from diffpy.Structure.atom import Atom
from diffpy.Structure.lattice import Lattice
from diffpy.Structure.structure import Structure
from diffpy.Structure.pdffitstructure import PDFFitStructure

# obtain version information
from diffpy.Structure.version import __version__

# top level routines

def loadStructure(filename, fmt='auto', **kw):
    """Load a new structure from a specified file.

    filename -- file to be loaded
    fmt      -- format of the structure file.  Must be one of the formats
                defined in the Parsers subpackage.  When 'auto', all
                Parsers are tried in sequence.
    kw       -- keyword arguments passed to the getParser factory function.

    Return a new Structure object.
    Return PDFFitStructure object for 'pdffit' or 'discus' formats.
    """
    from diffpy.Structure.Parsers import getParser
    p = getParser(fmt, **kw)
    rv = p.parseFile(filename)
    return rv

# silence pyflakes checker
assert (StructureFormatError and
        LatticeError and SymmetryError)
assert Atom
assert Lattice
assert Structure
assert PDFFitStructure
assert __version__ or True

# End of file

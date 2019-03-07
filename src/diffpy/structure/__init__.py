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

# Interface definitions ------------------------------------------------------

from diffpy.structure.structureerrors import (StructureFormatError,
        LatticeError, SymmetryError)
from diffpy.structure.atom import Atom
from diffpy.structure.lattice import Lattice
from diffpy.structure.structure import Structure
from diffpy.structure.pdffitstructure import PDFFitStructure

# obtain version information
from diffpy.structure.version import __version__

# top level routines

def loadStructure(filename, fmt='auto', **kw):
    """
    Load new structure object from the specified file.

    Parameters
    ----------

    filename : str
        Path to the file to be loaded.
    fmt : str, optional
        Format of the structure file such as 'cif' or 'xyz'.  Must be
        one of the formats listed by the `parsers.inputFormats` function.
        When 'auto', all supported formats are tried in a sequence.
    kw : misc, optional
        Extra keyword arguments that are passed to `parsers.getParser`
        function.  These configure the dedicated Parser object that
        is used to read content in filename.

    Returns
    -------
    stru : `Structure`, `PDFFitStructure`
        The new Structure object loaded from the specified file.
        Return a more specific PDFFitStructure type for 'pdffit'
        and 'discus' formats.
    """
    from diffpy.structure.parsers import getParser
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

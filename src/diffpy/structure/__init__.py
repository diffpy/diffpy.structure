#!/usr/bin/env python
##############################################################################
#
# (c) 2024-2026 The Trustees of Columbia University in the City of New York.
# All rights reserved.
#
# File coded by: Chris Farrow, Pavol Juhas,  Simon Billinge, Billinge Group members.
#
# See GitHub contributions for a more detailed list of contributors.
# https://github.com/diffpy/diffpy.structure/graphs/contributors
#
# See LICENSE.rst for license information.
#
##############################################################################
"""Crystal structure container and parsers for structure formats.

Classes related to the structure of materials:
    * Atom
    * Lattice
    * Structure
    * PDFFitStructure

Other classes:
    * SpaceGroup
    * SymOp
    * ExpandAsymmetricUnit
    * GeneratorSite
    * SymmetryConstraints

Exceptions:
    * StructureFormatError
    * LatticeError
    * SymmetryError
"""


import sys

import diffpy.structure as _structure
from diffpy.structure.atom import Atom
from diffpy.structure.lattice import Lattice
from diffpy.structure.parsers import get_parser
from diffpy.structure.pdffitstructure import PDFFitStructure
from diffpy.structure.structure import Structure
from diffpy.structure.structureerrors import LatticeError, StructureFormatError, SymmetryError

# package version
from diffpy.structure.version import __version__

# Deprecations -------------------------------------------------------


# @deprecated
# custom deprecator for diffpy.Structure module
class DeprecatedStructureModule:
    """Proxy for backward compatibility of diffpy.Structure."""

    def __getattr__(self, name):
        import warnings

        warnings.warn(
            "Module 'diffpy.Structure' is deprecated and will be removed in version 4.0. "
            "Use 'diffpy.structure' instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return getattr(_structure, name)


sys.modules["diffpy.Structure"] = DeprecatedStructureModule()

# top level routines


def loadStructure(filename, fmt="auto", **kw):
    """Load new structure object from the specified file.

    Parameters
    ----------

    filename : str
        Path to the file to be loaded.
    fmt : str, Optional
        Format of the structure file such as 'cif' or 'xyz'. Must be
        one of the formats listed by the `parsers.inputFormats` function.
        When 'auto', all supported formats are tried in a sequence.
    kw : Optional
        Extra keyword arguments that are passed to `parsers.getParser`
        function. These configure the dedicated Parser object that
        is used to read content in filename.

    Returns
    -------
    stru : `Structure`, `PDFFitStructure`
        The new Structure object loaded from the specified file.
        Return a more specific PDFFitStructure type for 'pdffit'
        and 'discus' formats.
    """

    p = get_parser(fmt, **kw)
    rv = p.parse_file(filename)
    return rv


# silence pyflakes checker
assert StructureFormatError and LatticeError and SymmetryError
assert Atom
assert Lattice
assert Structure
assert PDFFitStructure

# silence the pyflakes syntax checker
assert __version__ or True

# End of file

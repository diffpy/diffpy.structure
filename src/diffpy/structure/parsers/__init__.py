#!/usr/bin/env python
##############################################################################
#
# diffpy.structure  by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2007 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################

"""Conversion plugins for various structure formats.

The recognized structure formats are defined by subclassing `StructureParser`,
by convention these classes are named `P_<format>.py`. The parser classes should
to override the `parseLines()` and `toLines()` methods of `StructureParser`.
Any structure parser needs to be registered in `parser_index` module.

For normal usage it should be sufficient to use the routines provided
in this module.

Content:
    * StructureParser: base class for a concrete Parser
    * parser_index: dictionary of known structure formats
    * getParser: factory for Parser at given format
    * inputFormats: list of available input formats
    * outputFormats: list of available output formats
"""

from diffpy.structure.parsers.parser_index_mod import parser_index
from diffpy.structure.parsers.structureparser import StructureParser
from diffpy.structure.structureerrors import StructureFormatError

# silence pyflakes checker
assert StructureParser


def getParser(format, **kw):
    """Return Parser instance for a given structure format.

    Parameters
    ----------
    format : str
        String with the format name, see `parser_index_mod`.
    **kw : dict
        Keyword arguments passed to the Parser init function.

    Returns
    -------
    Parser
        Parser instance for the given format.

    Raises
    ------
    StructureFormatError
        When the format is not defined.
    """
    if format not in parser_index:
        emsg = "no parser for '%s' format" % format
        raise StructureFormatError(emsg)
    pmod = parser_index[format]["module"]
    ns = {}
    import_cmd = "from diffpy.structure.parsers import %s as pm" % pmod
    exec(import_cmd, ns)
    return ns["pm"].getParser(**kw)


def inputFormats():
    """Return list of implemented input structure formats."""
    input_formats = [fmt for fmt, prop in parser_index.items() if prop["has_input"]]
    input_formats.sort()
    return input_formats


def outputFormats():
    """Return list of implemented output structure formats."""
    output_formats = [fmt for fmt, prop in parser_index.items() if prop["has_output"]]
    output_formats.sort()
    return output_formats

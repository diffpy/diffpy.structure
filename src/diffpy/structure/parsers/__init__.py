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
from diffpy.utils._deprecator import build_deprecation_message, deprecated

# silence pyflakes checker
assert StructureParser

parsers_base = "diffpy.structure"
removal_version = "4.0.0"
getParser_deprecation_msg = build_deprecation_message(
    parsers_base,
    "getParser",
    "get_parser",
    removal_version,
)
inputFormats_deprecation_msg = build_deprecation_message(
    parsers_base,
    "inputFormats",
    "input_formats",
    removal_version,
)
outputFormats_deprecation_msg = build_deprecation_message(
    parsers_base,
    "outputFormats",
    "output_formats",
    removal_version,
)


@deprecated(getParser_deprecation_msg)
def getParser(format, **kw):
    """This function has been deprecated and will be removed in version
    4.0.0.

    Please use diffpy.structure.get_parser instead.
    """
    return get_parser(format, **kw)


def get_parser(format, **kw):
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
    return ns["pm"].get_parser(**kw)


@deprecated(inputFormats_deprecation_msg)
def inputFormats():
    """This function has been deprecated and will be removed in version
    4.0.0.

    Please use diffpy.structure.input_formats instead.
    """
    return input_formats()


def input_formats():
    """Return list of implemented input structure formats."""
    input_formats = [fmt for fmt, prop in parser_index.items() if prop["has_input"]]
    input_formats.sort()
    return input_formats


@deprecated(outputFormats_deprecation_msg)
def outputFormats():
    """This function has been deprecated and will be removed in version
    4.0.0.

    Please use diffpy.structure.output_formats instead.
    """
    return output_formats()


def output_formats():
    """Return list of implemented output structure formats."""
    output_formats = [fmt for fmt, prop in parser_index.items() if prop["has_output"]]
    output_formats.sort()
    return output_formats

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

"""Parser for automatic file format detection.

This Parser does not provide the the `toLines()` method.
"""

import os

from diffpy.structure.parsers import StructureParser, parser_index
from diffpy.structure.structureerrors import StructureFormatError


class P_auto(StructureParser):
    """Parser with automatic detection of structure format.

    This parser attempts to automatically detect the format of a given
    structure file and parse it accordingly. When successful, it sets
    its `format` attribute to the detected structure format.

    Parameters
    ----------
    **kw : dict
        Keyword arguments for the structure parser.

    Attributes
    ----------
    format : str
        Detected structure format. Initially set to "auto" and updated
        after successful detection of the structure format.
    pkw : dict
        Keyword arguments passed to the parser.
    """

    def __init__(self, **kw):
        StructureParser.__init__(self)
        self.format = "auto"
        self.pkw = kw
        return

    # parseLines helpers
    def _getOrderedFormats(self):
        """Build a list of relevance ordered structure formats.
        This only works when `self.filename` has a known extension.
        """
        from diffpy.structure.parsers import inputFormats

        ofmts = [fmt for fmt in inputFormats() if fmt != "auto"]
        if not self.filename:
            return ofmts
        # filename is defined here
        filebase = os.path.basename(self.filename)
        from fnmatch import fnmatch

        # loop over copy of ofmts
        for fmt in list(ofmts):
            pattern = parser_index[fmt]["file_pattern"]
            if pattern in ("*.*", "*"):
                continue
            anymatch = [1 for p in pattern.split("|") if fnmatch(filebase, p)]
            if anymatch:
                ofmts.remove(fmt)
                ofmts.insert(0, fmt)
        return ofmts

    def parseLines(self, lines):
        """Detect format and create `Structure` instance from a list of lines.

        Set format attribute to the detected file format.

        Parameters
        ----------
        lines : list
            List of lines with structure data.

        Returns
        -------
        Structure
            `Structure` object.

        Raises
        ------
        StructureFormatError
        """
        return self._wrapParseMethod("parseLines", lines)

    def parse(self, s):
        """Detect format and create `Structure` instance from a string.

        Set format attribute to the detected file format.

        Parameters
        ----------
        s : str
            String with structure data.

        Returns
        -------
        Structure
            `Structure` object.

        Raises
        ------
        StructureFormatError
        """
        return self._wrapParseMethod("parse", s)

    def parseFile(self, filename):
        """Detect format and create Structure instance from an existing file.

        Set format attribute to the detected file format.

        Parameters
        ----------
        filename : str
            Path to structure file.

        Returns
        -------
        Structure
            `Structure` object.

        Raises
        ------
        StructureFormatError
            If the structure format is unknown or invalid.
        IOError
            If the file cannot be read.
        """
        self.filename = filename
        return self._wrapParseMethod("parseFile", filename)

    def _wrapParseMethod(self, method, *args, **kwargs):
        """A helper evaluator method that try the specified parse method with
        each registered structure parser and return the first successful
        resul.

        Structure parsers that match structure file extension are
        tried first.

        Parameters
        ----------
        method : str
            Name of the parse method to call.
        *args : tuple
            Positional arguments for the parse method.
        **kwargs : dict
            Keyword arguments for the parse method.

        Returns
        -------
        Structure
            `Structure` object.

        Raises
        ------
        StructureFormatError
        """
        from diffpy.structure.parsers import getParser

        ofmts = self._getOrderedFormats()
        stru = None
        # try all parsers in sequence
        parsers_emsgs = []
        for fmt in ofmts:
            p = getParser(fmt, **self.pkw)
            try:
                pmethod = getattr(p, method)
                stru = pmethod(*args, **kwargs)
                self.format = fmt
                break
            except StructureFormatError as err:
                parsers_emsgs.append("%s: %s" % (fmt, err))
            except NotImplementedError:
                pass
        if stru is None:
            emsg = "\n".join(
                ["Unknown or invalid structure format.", "Errors per each tested structure format:"]
                + parsers_emsgs
            )
            raise StructureFormatError(emsg)
        self.__dict__.update(p.__dict__)
        return stru


# End of class P_auto

# Routines -------------------------------------------------------------------


def getParser(**kw):
    """Return a new instance of the automatic parser.

    Parameters
    ----------
    **kw : dict
        Keyword arguments for the structure parser

    Returns
    -------
    P_auto
        Instance of `P_auto`.
    """
    return P_auto(**kw)

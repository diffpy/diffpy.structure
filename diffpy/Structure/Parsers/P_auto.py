#!/usr/bin/env python
##############################################################################
#
# diffpy.Structure  by DANSE Diffraction group
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

"""Parser for automatic file format detection

This Parser does not provide the the toLines() method.
"""

import os

from diffpy.Structure import StructureFormatError
from diffpy.Structure.Parsers import StructureParser
from diffpy.Structure.Parsers import parser_index

class P_auto(StructureParser):
    """Parser with automatic detection of structure format.
    When successful, it sets its format attribute to detected
    structure format.
    """

    def __init__(self, **kw):
        StructureParser.__init__(self)
        self.format = "auto"
        self.pkw = kw
        return

    # parseLines helpers
    def _getOrderedFormats(self):
        """Build a list of relevance ordered structure formats.
        This only works when self.filename has a known extension.
        """
        from diffpy.Structure.Parsers import inputFormats
        ofmts = [fmt for fmt in inputFormats() if fmt != 'auto']
        if not self.filename:   return ofmts
        # filename is defined here
        filebase = os.path.basename(self.filename)
        from fnmatch import fnmatch
        # loop over copy of ofmts
        for fmt in list(ofmts):
            pattern = parser_index[fmt]['file_pattern']
            if pattern in ('*.*', '*'):     continue
            anymatch = [1 for p in pattern.split('|') if fnmatch(filebase, p)]
            if anymatch:
                ofmts.remove(fmt)
                ofmts.insert(0, fmt)
        return ofmts


    def parseLines(self, lines):
        """Detect format and create Structure instance from a list of lines.
        Set format attribute to the detected file format.

        Return Structure object or raise StructureFormatError exception.
        """
        return self._wrapParseMethod("parseLines", lines)


    def parse(self, s):
        """Detect format and create Structure instance from a string.
        Set format attribute to the detected file format.

        Return Structure object or raise StructureFormatError exception.
        """
        return self._wrapParseMethod('parse', s)


    def parseFile(self, filename):
        '''Detect format and create Structure instance from an existing file.
        Set format attribute to the detected file format.

        filename  -- path to structure file

        Return Structure object.
        Raise StructureFormatError or IOError.
        '''
        self.filename = filename
        return self._wrapParseMethod("parseFile", filename)


    def _wrapParseMethod(self, method, *args, **kwargs):
        """A helper evaluator method.  Try the specified parse method with
        each registered structure parser and return the first successful
        resul.  Structure parsers that match structure file extension are
        tried first.

        Set format attribute to the detected file format.
        Return Structure instance, or raise StructureFormatError.
        """
        from diffpy.Structure.Parsers import getParser
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
            except StructureFormatError, err:
                parsers_emsgs.append("%s: %s" % (fmt, err))
            except NotImplementedError:
                pass
        if stru is None:
            emsg = "\n".join([
                "Unknown or invalid structure format.",
                "Errors per each tested structure format:"] + parsers_emsgs)
            raise StructureFormatError(emsg)
        self.__dict__.update(p.__dict__)
        return stru
    # End of parseLines

# End of class P_auto

# Routines

def getParser(**kw):
    return P_auto(**kw)

# End of file

########################################################################
#
# Structure         by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2007 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
########################################################################

"""Parser for automatic file format detection

This Parser does not provide the the toLines() method.
"""

__id__ = "$Id$"

from import_helper import InvalidStructureFormat
from StructureParser import StructureParser
from parser_index import parser_index

class P_auto(StructureParser):
    """Parser with automatic detection of structure format.
    When successful, it sets its format attribute to detected
    structure format.
    """

    def __init__(self):
        StructureParser.__init__(self)
        self.format = "auto"
        return

    # parseLines helpers
    def getOrderedFormats(self):
        """Build a list of relevance ordered structure formats.
        This only works when self.filename has a known extension.
        """
        import os.path
        from __init__ import inputFormats
        ofmts = [fmt in inputFormats() if fmt != 'auto']
        if not self.filename:   return ofmts
        # filename is defined here
        filebase = os.path.basename(self.filename)
        from fnmatch import fnmatch
        # loop over copy of ofmts
        for fmt in list(ofmts):
            pattern = parser_index[fmt]['file_pattern']
            if pattern in ('*.*', '*'):     continue
            anymatch = [1 for p in pattern.split('|') if fnmatch(basename, p)]
            if anymatch:
                ofmts.remove(fmt)
                ofmts.insert(0, fmt)
        return ofmts

    def parseLines(self, lines):
        """Detect format and parse given list of lines.
        Set format attribute to the detected file format.

        Return Structure instance, or raise InvalidStructureFormat.
        """
        ofmts = self.getOrderedFormats()
        from __init__ import getParser
        stru = None
        # try all parsers in sequence
        for fmt in ofmts:
            p = getParser(fmt)
            try:
                stru = p.parseLines(lines)
                self.format = fmt
                break
            except (InvalidStructureFormat, NotImplementedError):
                pass
        if stru is None:
            raise InvalidStructureFormat, "unknown structure format"
        return stru
    # End of parseLines

# End of class P_auto

# Routines

def getParser():
    return P_auto()

# End of file

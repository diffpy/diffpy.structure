#!/usr/bin/env python
##############################################################################
#
# diffpy.structure  by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2008 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################

"""Definition of StructureParser, a base class for specific parsers.
"""


class StructureParser(object):
    """Base class for all structure parsers.

    Attributes
    ----------
    format : str
        Format name of particular parser.
    filename : str
        Path to structure file that is read or written.
    """

    def __init__(self):
        self.format = None
        self.filename = None
        return

    def parseLines(self, lines):
        """Create Structure instance from a list of lines.

        Return Structure object or raise StructureFormatError exception.

        Note
        ----
        This method has to be overloaded in derived class.
        """
        raise NotImplementedError("parseLines not defined for '%s' format" % self.format)
        return

    def toLines(self, stru):
        """Convert Structure stru to a list of lines.

        Return list of strings.

        Note
        ----
        This method has to be overloaded in derived class.
        """
        raise NotImplementedError("toLines not defined for '%s' format" % self.format)

    def parse(self, s):
        """Create `Structure` instance from a string."""
        lines = s.rstrip("\r\n").split("\n")
        stru = self.parseLines(lines)
        return stru

    def tostring(self, stru):
        """Convert `Structure` instance to a string."""
        lines = self.toLines(stru)
        s = "\n".join(lines) + "\n"
        return s

    def parseFile(self, filename):
        """Create Structure instance from an existing file."""
        self.filename = filename
        with open(filename) as fp:
            s = fp.read()
        stru = self.parse(s)
        return stru


# End of class StructureParser

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

"""definition of StructureParser superclass"""

__id__ = "$Id$"

from Structure.exceptions import InvalidStructureFormat

##############################################################################
class StructureParser:
    """StructureParser --> superclass for structure parsers

    Data members:
        format -- format of particular parser
    """

    def __init__(self):
        self.format = None
        return

    def parseLines(self, lines):
        """parse list of lines obtained from structure file

        return Structure object or raise InvalidStructureFormat exception.
        This method has to be overloaded in a derived class
        """
        raise NotImplementedError, \
                "parseLines not defined for '%s' format" % self.format
        return

    def toLines(self, stru):
        """Convert Structure stru to a list of lines.
        This method has to be overloaded in a derived class.

        Return list of strings.
        """
        raise NotImplementedError, \
                "toLines not defined for '%s' format" % self.format

# End of StructureParser

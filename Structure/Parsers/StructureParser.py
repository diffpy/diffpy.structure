"""definition of StructureParser superclass"""

__id__ = "$Id$"

from Structure.structure import InvalidStructureFormat

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
        """convert Structure stru to a list of lines
        
        return a list of strings. 
        This method has to be overloaded in a derived class.
        """
        raise NotImplementedError, \
                "toLines not defined for '%s' format" % self.format

# End of StructureParser

##############################################################################
# helper function

def isfloat(s):
    """True if argument can be converted to float"""
    try:
        x = float(s)
        return True
    except ValueError:
        return False

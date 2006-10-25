"""Parser with automatic file format detection

This Parser leaves undefined the toLines() method
"""

__id__ = "$Id$"

from Structure.exceptions import InvalidStructureFormat
from StructureParser import StructureParser

class Parser(StructureParser):
    """Parser --> StructureParser subclass for automatic format"""

    def __init__(self):
        self.format = "auto"
        self.lastformat = None
        return

    def parseLines(self, lines):
        """Detect format and parse given list of lines.

        Return Structure instance, or raise InvalidStructureFormat.
        """
        stru = None
        # try all parsers in sequence
        import Structure.Parsers as Parsers
        mods = {}
        mods.update(Parsers.modules)
        del mods[self.format]
        for pfmt, pmod in mods.iteritems():
            try:
                exec('from ' + pmod + ' import Parser')
                p = Parser()
                stru = p.parseLines(lines)
                self.lastformat = pfmt
                break
            except (InvalidStructureFormat, NotImplementedError):
                pass
        if stru is None:
            raise InvalidStructureFormat, "unknown structure format"
        return stru
    # End of parseLines

# End of Parser

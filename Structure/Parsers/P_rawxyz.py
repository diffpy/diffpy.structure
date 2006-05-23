"""Parser for raw XYZ file format"""

__id__ = "$Id$"

import sys
from Structure.structure import Structure, InvalidStructureFormat
from Structure.lattice import Lattice
from Structure.atom import Atom
from StructureParser import StructureParser, isfloat

class Parser(StructureParser):
    """Parser --> StructureParser subclass for RAWXYZ format"""

    def __init__(self):
        self.format = "rawxyz"
        return

    def parseLines(self, lines):
        """parse list of lines in RAWXYZ format

        return Structure object or raise InvalidStructureFormat exception
        """
        linefields = [l.split() for l in lines]
        # prepare output structure
        stru = Structure()
        # find first valid record
        start = 0
        for field in linefields:
            if len(field) == 0 or field[0] == "#":
                start += 1
            else:
                break
        # find the last valid record
        stop = len(lines)
        while stop > start and len(linefields[stop-1]) == 0:
            stop -= 1
        # get out for empty structure
        if start >= stop:
            return stru
        # here we have at least one valid record line
        # figure out xyz layout from the first line for plain and raw formats
        floatfields = [ isfloat(f) for f in linefields[start] ]
        nfields = len(linefields[start])
        if nfields not in (3, 4):
            raise InvalidStructureFormat, ("%d: invalid RAWXYZ format, " +
                    "expected 3 or 4 columns") % (start+1)
        if floatfields[:3] == [True, True, True]:
            el_idx, x_idx = (None, 0)
        elif floatfields[:4] == [False, True, True, True]:
            el_idx, x_idx = (0, 1)
        else:
            raise InvalidStructureFormat, \
                    "%d: invalid RAWXYZ format" % (start+1)
        # now try to read all record lines
        try:
            p_nl = start
            for fields in linefields[start:] :
                p_nl += 1
                if fields == []:
                    continue
                elif len(fields) != nfields:
                    raise InvalidStructureFormat, ('%d: all lines must have ' +
                            'the same number of columns') % p_nl
                element = el_idx is not None and fields[el_idx] or ""
                xyz = [ float(f) for f in fields[x_idx:x_idx+3] ]
                if len(xyz) == 2:
                    xyz.append(0.0)
                stru.append(Atom(element, xyz=xyz))
        except ValueError:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            raise InvalidStructureFormat, \
                    "%d: invalid number" % p_nl, exc_traceback
        return stru
    # End of parseLines

    def toLines(self, stru):
        """convert Structure stru to a list of lines in XYZ format"""
        lines = []
        for a in stru:
            rc = stru.cartesian(a)
            s = "%s %g %g %g" % (a.element, rc[0], rc[1], rc[2])
            lines.append(s.lstrip())
        return lines
    # End of toLines

# End of Parser

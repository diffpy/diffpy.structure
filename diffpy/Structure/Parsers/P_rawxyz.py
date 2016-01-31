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

"""Parser for raw XYZ file format.  Raw XYZ is a 3 or 4 column text
file with cartesian coordinates of atoms and an optional first column
for atom types.
"""

import sys

from diffpy.Structure import Structure
from diffpy.Structure import StructureFormatError
from diffpy.Structure.utils import isfloat
from diffpy.Structure.Parsers import StructureParser

class P_rawxyz(StructureParser):
    """Parser --> StructureParser subclass for RAWXYZ format"""

    def __init__(self):
        StructureParser.__init__(self)
        self.format = "rawxyz"
        return

    def parseLines(self, lines):
        """Parse list of lines in RAWXYZ format.

        Return Structure object or raise StructureFormatError.
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
            emsg = ("%d: invalid RAWXYZ format, expected 3 or 4 columns" %
                    (start + 1))
            raise StructureFormatError(emsg)
        if floatfields[:3] == [True, True, True]:
            el_idx, x_idx = (None, 0)
        elif floatfields[:4] == [False, True, True, True]:
            el_idx, x_idx = (0, 1)
        else:
            emsg = "%d: invalid RAWXYZ format" % (start + 1)
            raise StructureFormatError(emsg)
        # now try to read all record lines
        try:
            p_nl = start
            for fields in linefields[start:] :
                p_nl += 1
                if fields == []:
                    continue
                elif len(fields) != nfields:
                    emsg = ('%d: all lines must have ' +
                            'the same number of columns') % p_nl
                    raise StructureFormatError, emsg
                element = el_idx is not None and fields[el_idx] or ""
                xyz = [ float(f) for f in fields[x_idx:x_idx+3] ]
                if len(xyz) == 2:
                    xyz.append(0.0)
                stru.addNewAtom(element, xyz=xyz)
        except ValueError:
            emsg = "%d: invalid number" % p_nl
            exc_type, exc_value, exc_traceback = sys.exc_info()
            raise StructureFormatError, emsg, exc_traceback
        return stru
    # End of parseLines

    def toLines(self, stru):
        """Convert Structure stru to a list of lines in XYZ format.

        Return list of strings.
        """
        lines = []
        for a in stru:
            rc = a.xyz_cartn
            s = "%s %g %g %g" % (a.element, rc[0], rc[1], rc[2])
            lines.append(s.lstrip())
        return lines
    # End of toLines

# End of class P_rawxyz

# Routines

def getParser():
    return P_rawxyz()

# End of file

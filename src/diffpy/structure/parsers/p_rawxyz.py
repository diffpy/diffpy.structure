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

"""Parser for raw XYZ file format.

Raw XYZ is a 3 or 4 column text file with cartesian coordinates
of atoms and an optional first column for atom types.
"""

import sys

from diffpy.structure import Structure
from diffpy.structure.parsers import StructureParser
from diffpy.structure.structureerrors import StructureFormatError
from diffpy.structure.utils import isfloat


class P_rawxyz(StructureParser):
    """Parser --> StructureParser subclass for RAWXYZ format.

    Attributes
    ----------
    format : str
        Format name, default "rawxyz".
    """

    def __init__(self):
        StructureParser.__init__(self)
        self.format = "rawxyz"
        return

    def parseLines(self, lines):
        """Parse list of lines in RAWXYZ format.

        Parameters
        ----------
        lines : list of str
            List of lines in RAWXYZ format.

        Returns
        -------
        Structure
            Parsed structure instance.

        Raises
        ------
        StructureFormatError
            Invalid RAWXYZ format.
        """
        linefields = [line.split() for line in lines]
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
        while stop > start and len(linefields[stop - 1]) == 0:
            stop -= 1
        # get out for empty structure
        if start >= stop:
            return stru
        # here we have at least one valid record line
        # figure out xyz layout from the first line for plain and raw formats
        floatfields = [isfloat(f) for f in linefields[start]]
        nfields = len(linefields[start])
        if nfields not in (3, 4):
            emsg = "%d: invalid RAWXYZ format, expected 3 or 4 columns" % (start + 1)
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
            for fields in linefields[start:]:
                p_nl += 1
                if fields == []:
                    continue
                elif len(fields) != nfields:
                    emsg = ("%d: all lines must have " + "the same number of columns") % p_nl
                    raise StructureFormatError(emsg)
                element = el_idx is not None and fields[el_idx] or ""
                xyz = [float(f) for f in fields[x_idx : x_idx + 3]]
                if len(xyz) == 2:
                    xyz.append(0.0)
                stru.addNewAtom(element, xyz=xyz)
        except ValueError:
            emsg = "%d: invalid number" % p_nl
            exc_type, exc_value, exc_traceback = sys.exc_info()
            e = StructureFormatError(emsg)
            raise e.with_traceback(exc_traceback)
        return stru

    def toLines(self, stru):
        """Convert Structure stru to a list of lines in RAWXYZ format.

        Parameters
        ----------
        stru : Structure
            Structure to be converted.

        Returns
        -------
        list of str
            List of lines in RAWXYZ format.
        """
        lines = []
        for a in stru:
            rc = a.xyz_cartn
            s = "%s %g %g %g" % (a.element, rc[0], rc[1], rc[2])
            lines.append(s.lstrip())
        return lines


# End of class P_rawxyz

# Routines -------------------------------------------------------------------


def getParser():
    """Return new `parser` object for RAWXYZ format.

    Returns
    -------
    P_rawxyz
        Instance of `P_rawxyz`.
    """
    return P_rawxyz()

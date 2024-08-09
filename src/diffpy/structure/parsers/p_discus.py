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

"""Parser for DISCUS structure format
"""

import sys
from functools import reduce

from diffpy.structure import Lattice, PDFFitStructure
from diffpy.structure.parsers import StructureParser
from diffpy.structure.structureerrors import StructureFormatError


class P_discus(StructureParser):
    """Parser for DISCUS structure format. The parser chokes
    on molecule and generator records.

    Attributes
    ----------
    format : str
        File format name, default "discus".
    nl : int
        Line number of the current line being parsed.
    lines : list of str
        List of lines from the input file.
    line : str
        Current line being parsed.
    stru : PDFFitStructure
        Structure being parsed.
    ignored_lines : list of str
        List of lines that were ignored during parsing.
    cell_read : bool
        ``True`` if cell record processed.
    ncell_read : bool
        ``True`` if ncell record processed.
    """

    def __init__(self):
        StructureParser.__init__(self)
        self.format = "discus"
        # helper variables
        self.nl = None
        self.lines = None
        self.line = None
        self.stru = None
        self.ignored_lines = []
        self.cell_read = False
        self.ncell_read = False
        return

    def parseLines(self, lines):
        """Parse list of lines in DISCUS format.

        Parameters
        ----------
        lines : list of str
            List of lines from the input file.

        Returns
        -------
        PDFFitStructure
            Parsed `PDFFitStructure` instance.

        Raises
        ------
        StructureFormatError
            If the file is not in DISCUS format.
        """
        self.lines = lines
        ilines = self._linesIterator()
        self.stru = PDFFitStructure()
        record_parsers = {
            "cell": self._parse_cell,
            "format": self._parse_format,
            "generator": self._parse_not_implemented,
            "molecule": self._parse_not_implemented,
            "ncell": self._parse_ncell,
            "spcgr": self._parse_spcgr,
            "symmetry": self._parse_not_implemented,
            "title": self._parse_title,
            "shape": self._parse_shape,
        }
        try:
            # parse header
            for self.line in ilines:
                words = self.line.split()
                if not words or words[0][0] == "#":
                    continue
                if words[0] == "atoms":
                    break
                rp = record_parsers.get(words[0], self._parse_unknown_record)
                rp(words)
            # check if cell has been defined
            if not self.cell_read:
                emsg = "%d: unit cell not defined" % self.nl
                raise StructureFormatError(emsg)
            # parse atoms
            for self.line in ilines:
                words = self.line.replace(",", " ").split()
                if not words or words[0][0] == "#":
                    continue
                self._parse_atom(words)
            # self consistency check
            exp_natoms = reduce(lambda x, y: x * y, self.stru.pdffit["ncell"])
            # only check if ncell record exists
            if self.ncell_read and exp_natoms != len(self.stru):
                emsg = "Expected %d atoms, read %d." % (exp_natoms, len(self.stru))
                raise StructureFormatError(emsg)
            # take care of superlattice
            if self.stru.pdffit["ncell"][:3] != [1, 1, 1]:
                latpars = list(self.stru.lattice.abcABG())
                superlatpars = [latpars[i] * self.stru.pdffit["ncell"][i] for i in range(3)] + latpars[3:]
                superlattice = Lattice(*superlatpars)
                self.stru.placeInLattice(superlattice)
                self.stru.pdffit["ncell"] = [1, 1, 1, exp_natoms]
        except (ValueError, IndexError):
            exc_type, exc_value, exc_traceback = sys.exc_info()
            emsg = "%d: file is not in DISCUS format" % self.nl
            e = StructureFormatError(emsg)
            raise e.with_traceback(exc_traceback)
        return self.stru

    def toLines(self, stru):
        """Convert `Structure` stru to a list of lines in DISCUS format.

        Parameters
        ----------
        stru : Structure
            Structure to be converted.

        Returns
        -------
        list of str
            List of lines in DISCUS format.
        """
        self.stru = stru
        # if necessary, convert self.stru to PDFFitStructure
        if not isinstance(stru, PDFFitStructure):
            self.stru = PDFFitStructure(stru)
        # build the stru_pdffit dictionary initialized from the defaults
        # in PDFFitStructure
        stru_pdffit = PDFFitStructure().pdffit
        if stru.pdffit:
            stru_pdffit.update(stru.pdffit)
        # here we can start
        self.lines = lines = []
        lines.append(("title   " + self.stru.title).strip())
        lines.append("spcgr   " + stru_pdffit["spcgr"])
        if stru_pdffit.get("spdiameter", 0.0) > 0.0:
            line = "shape   sphere, %g" % stru_pdffit["spdiameter"]
            lines.append(line)
        if stru_pdffit.get("stepcut", 0.0) > 0.0:
            line = "shape   stepcut, %g" % stru_pdffit["stepcut"]
            lines.append(line)
        lines.append("cell   %9.6f, %9.6f, %9.6f, %9.6f, %9.6f, %9.6f" % self.stru.lattice.abcABG())
        lines.append("ncell  %9i, %9i, %9i, %9i" % (1, 1, 1, len(self.stru)))
        lines.append("atoms")
        for a in self.stru:
            lines.append(
                "%-4s %17.8f %17.8f %17.8f %12.4f" % (a.element.upper(), a.xyz[0], a.xyz[1], a.xyz[2], a.Bisoequiv)
            )
        return lines

    def _linesIterator(self):
        """Iterator over `self.lines`, which increments `self.nl`"""
        # ignore trailing empty lines
        stop = len(self.lines)
        while stop > 0 and self.lines[stop - 1].strip() == "":
            stop -= 1
        self.nl = 0
        # read header of PDFFit file
        for self.line in self.lines[:stop]:
            self.nl += 1
            yield self.line
        pass

    def _parse_cell(self, words):
        """Process the cell record from DISCUS structure file."""
        # split again on spaces or commas
        words = self.line.replace(",", " ").split()
        latpars = [float(w) for w in words[1:7]]
        try:
            self.stru.lattice.setLatPar(*latpars)
        except ZeroDivisionError:
            emsg = "%d: Invalid lattice parameters - zero cell volume" % self.nl
            raise StructureFormatError(emsg)
        self.cell_read = True
        return

    def _parse_format(self, words):
        """Process the format record from DISCUS structure file."""
        if words[1] == "pdffit":
            emsg = "%d: file is not in DISCUS format" % self.nl
            raise StructureFormatError(emsg)
        return

    def _parse_ncell(self, words):
        """Process the ncell record from DISCUS structure file."""
        # split again on spaces or commas
        words = self.line.replace(",", " ").split()
        self.stru.pdffit["ncell"] = [int(w) for w in words[1:5]]
        self.ncell_read = True
        return

    def _parse_spcgr(self, words):
        """Process the spcgr record from DISCUS structure file."""
        self.stru.pdffit["spcgr"] = "".join(words[1:])
        return

    def _parse_title(self, words):
        """Process the title record from DISCUS structure file."""
        self.stru.title = self.line.lstrip()[5:].strip()
        return

    def _parse_shape(self, words):
        """Process the shape record from DISCUS structure file.

        Parameters
        ----------
        words : list of str
            List of words in the line.

        Raises
        ------
        StructureFormatError
            Invalid type of particle shape correction.
        """
        # strip away any commas
        linefixed = " ".join(words).replace(",", " ")
        wordsfixed = linefixed.split()
        shapetype = wordsfixed[1]
        if shapetype == "sphere":
            self.stru.pdffit["spdiameter"] = float(words[2])
        elif shapetype == "stepcut":
            self.stru.pdffit["stepcut"] = float(words[2])
        else:
            emsg = "Invalid type of particle shape correction %r" % shapetype
            raise StructureFormatError(emsg)
        return

    def _parse_atom(self, words):
        """Process atom records in DISCUS structure file."""
        element = words[0][0:1].upper() + words[0][1:].lower()
        xyz = [float(w) for w in words[1:4]]
        Biso = float(words[4])
        self.stru.addNewAtom(element, xyz)
        a = self.stru.getLastAtom()
        a.Bisoequiv = Biso
        return

    def _parse_unknown_record(self, words):
        """Process unknown record in DISCUS structure file.

        Silently ignores the line and adds it to `self.ignored_lines`.

        Parameters
        ----------
        words : list of str
            List of words in the line.

        Raises
        ------
        StructureFormatError
            Unkown record.
        """
        self.ignored_lines.append(self.line)
        return

    def _parse_not_implemented(self, words):
        """Process the unimplemented records from DISCUS structure file.

        Parameters
        ----------
        words : list of str
            List of words in the line.

        Raises
        ------
        NotImplementedError
            If the record is not implemented.
        """
        emsg = "%d: reading of DISCUS record %r is not implemented." % (self.nl, words[0])
        raise NotImplementedError(emsg)


# End of class P_pdffit

# Routines -------------------------------------------------------------------


def getParser():
    """Return new `parser` object for DISCUS format.

    Returns
    -------
    P_discus
        Instance of `P_discus`.
    """
    return P_discus()

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

"""Parser for PDFfit structure format
"""

import sys
from functools import reduce

import numpy

from diffpy.structure import Lattice, PDFFitStructure
from diffpy.structure.parsers import StructureParser
from diffpy.structure.structureerrors import StructureFormatError


class P_pdffit(StructureParser):
    """Parser for PDFfit structure format.

    Attributes
    ----------
    format : str
        Format name, default "pdffit".
    ignored_lines : list
        List of lines ignored during parsing.
    stru : PDFFitStructure
        Structure instance used for cif input or output.
    """

    def __init__(self):
        StructureParser.__init__(self)
        self.format = "pdffit"
        self.ignored_lines = []
        self.stru = None
        return

    def parseLines(self, lines):
        """Parse list of lines in PDFfit format.

        Parameters
        ----------
        lines : list of str
            List of lines in PDB format.

        Returns
        -------
        Structure
            Parsed structure instance.

        Raises
        ------
        StructureFormatError
            File not in PDFfit format.
        """
        p_nl = 0
        try:
            self.stru = PDFFitStructure()
            stru = self.stru
            cell_line_read = False
            stop = len(lines)
            while stop > 0 and lines[stop - 1].strip() == "":
                stop -= 1
            ilines = iter(lines[:stop])
            # read header of PDFFit file
            for line in ilines:
                p_nl += 1
                words = line.split()
                if len(words) == 0 or words[0][0] == "#":
                    continue
                elif words[0] == "title":
                    stru.title = line.lstrip()[5:].strip()
                elif words[0] == "scale":
                    stru.pdffit["scale"] = float(words[1])
                elif words[0] == "sharp":
                    l1 = line.replace(",", " ")
                    sharp_pars = [float(w) for w in l1.split()[1:]]
                    if len(sharp_pars) < 4:
                        stru.pdffit["delta2"] = sharp_pars[0]
                        stru.pdffit["sratio"] = sharp_pars[1]
                        stru.pdffit["rcut"] = sharp_pars[2]
                    else:
                        stru.pdffit["delta2"] = sharp_pars[0]
                        stru.pdffit["delta1"] = sharp_pars[1]
                        stru.pdffit["sratio"] = sharp_pars[2]
                        stru.pdffit["rcut"] = sharp_pars[3]
                elif words[0] == "spcgr":
                    key = "spcgr"
                    start = line.find(key) + len(key)
                    value = line[start:].strip()
                    stru.pdffit["spcgr"] = value
                elif words[0] == "shape":
                    self._parse_shape(line)
                elif words[0] == "cell":
                    cell_line_read = True
                    l1 = line.replace(",", " ")
                    latpars = [float(w) for w in l1.split()[1:7]]
                    stru.lattice = Lattice(*latpars)
                elif words[0] == "dcell":
                    l1 = line.replace(",", " ")
                    stru.pdffit["dcell"] = [float(w) for w in l1.split()[1:7]]
                elif words[0] == "ncell":
                    l1 = line.replace(",", " ")
                    stru.pdffit["ncell"] = [int(w) for w in l1.split()[1:5]]
                elif words[0] == "format":
                    if words[1] != "pdffit":
                        emsg = "%d: file is not in PDFfit format" % p_nl
                        raise StructureFormatError(emsg)
                elif words[0] == "atoms" and cell_line_read:
                    break
                else:
                    self.ignored_lines.append(line)
            # Header reading finished, check if required lines were present.
            if not cell_line_read:
                emsg = "%d: file is not in PDFfit format" % p_nl
                raise StructureFormatError(emsg)
            # Load data from atom entries.
            p_natoms = reduce(lambda x, y: x * y, stru.pdffit["ncell"])
            # we are now inside data block
            for line in ilines:
                p_nl += 1
                wl1 = line.split()
                element = wl1[0][0].upper() + wl1[0][1:].lower()
                xyz = [float(w) for w in wl1[1:4]]
                occ = float(wl1[4])
                stru.addNewAtom(element, xyz=xyz, occupancy=occ)
                a = stru.getLastAtom()
                p_nl += 1
                wl2 = next(ilines).split()
                a.sigxyz = [float(w) for w in wl2[0:3]]
                a.sigo = float(wl2[3])
                p_nl += 1
                wl3 = next(ilines).split()
                p_nl += 1
                wl4 = next(ilines).split()
                p_nl += 1
                wl5 = next(ilines).split()
                p_nl += 1
                wl6 = next(ilines).split()
                U = numpy.zeros((3, 3), dtype=float)
                sigU = numpy.zeros((3, 3), dtype=float)
                U[0, 0] = float(wl3[0])
                U[1, 1] = float(wl3[1])
                U[2, 2] = float(wl3[2])
                sigU[0, 0] = float(wl4[0])
                sigU[1, 1] = float(wl4[1])
                sigU[2, 2] = float(wl4[2])
                U[0, 1] = U[1, 0] = float(wl5[0])
                U[0, 2] = U[2, 0] = float(wl5[1])
                U[1, 2] = U[2, 1] = float(wl5[2])
                sigU[0, 1] = sigU[1, 0] = float(wl6[0])
                sigU[0, 2] = sigU[2, 0] = float(wl6[1])
                sigU[1, 2] = sigU[2, 1] = float(wl6[2])
                a.anisotropy = stru.lattice.isanisotropic(U)
                a.U = U
                a.sigU = sigU
            if len(stru) != p_natoms:
                emsg = "expected %d atoms, read %d" % (p_natoms, len(stru))
                raise StructureFormatError(emsg)
            if stru.pdffit["ncell"][:3] != [1, 1, 1]:
                superlatpars = [latpars[i] * stru.pdffit["ncell"][i] for i in range(3)] + latpars[3:]
                superlattice = Lattice(*superlatpars)
                stru.placeInLattice(superlattice)
                stru.pdffit["ncell"] = [1, 1, 1, p_natoms]
        except (ValueError, IndexError):
            emsg = "%d: file is not in PDFfit format" % p_nl
            exc_type, exc_value, exc_traceback = sys.exc_info()
            e = StructureFormatError(emsg)
            raise e.with_traceback(exc_traceback)
        return stru

    def toLines(self, stru):
        """Convert `Structure` stru to a list of lines in PDFfit format.

        Parameters
        ----------
        stru : Structure
            Structure to be converted.

        Returns
        -------
        list of str
            List of lines in PDFfit format.
        """
        # build the stru_pdffit dictionary initialized from the defaults
        # in PDFFitStructure
        stru_pdffit = PDFFitStructure().pdffit
        if stru.pdffit:
            stru_pdffit.update(stru.pdffit)
        lines = []
        # default values of standard deviations
        d_sigxyz = numpy.zeros(3, dtype=float)
        d_sigo = 0.0
        d_sigU = numpy.zeros((3, 3), dtype=float)
        # here we can start
        line = "title  " + stru.title
        lines.append(line.strip())
        lines.append("format pdffit")
        lines.append("scale  %9.6f" % stru_pdffit["scale"])
        lines.append(
            "sharp  %9.6f, %9.6f, %9.6f, %9.6f"
            % (stru_pdffit["delta2"], stru_pdffit["delta1"], stru_pdffit["sratio"], stru_pdffit["rcut"])
        )
        lines.append("spcgr   " + stru_pdffit["spcgr"])
        if stru_pdffit.get("spdiameter", 0.0) > 0.0:
            line = "shape   sphere, %g" % stru_pdffit["spdiameter"]
            lines.append(line)
        if stru_pdffit.get("stepcut", 0.0) > 0.0:
            line = "shape   stepcut, %g" % stru_pdffit["stepcut"]
            lines.append(line)
        lat = stru.lattice
        lines.append(
            "cell   %9.6f, %9.6f, %9.6f, %9.6f, %9.6f, %9.6f"
            % (lat.a, lat.b, lat.c, lat.alpha, lat.beta, lat.gamma)
        )
        lines.append("dcell  %9.6f, %9.6f, %9.6f, %9.6f, %9.6f, %9.6f" % tuple(stru_pdffit["dcell"]))
        lines.append("ncell  %9i, %9i, %9i, %9i" % (1, 1, 1, len(stru)))
        lines.append("atoms")
        for a in stru:
            ad = a.__dict__
            lines.append(
                "%-4s %17.8f %17.8f %17.8f %12.4f" % (a.element.upper(), a.xyz[0], a.xyz[1], a.xyz[2], a.occupancy)
            )
            sigmas = numpy.concatenate((ad.get("sigxyz", d_sigxyz), [ad.get("sigo", d_sigo)]))
            lines.append("    %18.8f %17.8f %17.8f %12.4f" % tuple(sigmas))
            sigU = ad.get("sigU", d_sigU)
            Uii = (a.U[0][0], a.U[1][1], a.U[2][2])
            Uij = (a.U[0][1], a.U[0][2], a.U[1][2])
            sigUii = (sigU[0][0], sigU[1][1], sigU[2][2])
            sigUij = (sigU[0][1], sigU[0][2], sigU[1][2])
            lines.append("    %18.8f %17.8f %17.8f" % Uii)
            lines.append("    %18.8f %17.8f %17.8f" % sigUii)
            lines.append("    %18.8f %17.8f %17.8f" % Uij)
            lines.append("    %18.8f %17.8f %17.8f" % sigUij)
        return lines

    # Protected methods ------------------------------------------------------

    def _parse_shape(self, line):
        """Process shape line from PDFfit file and update self.stru.

        Parameters
        ----------
        line : str
            Line containing data for particle shape correction.

        Raises
        ------
        StructureFormatError
            Invalid type of particle shape correction.
        """
        line_nocommas = line.replace(",", " ")
        words = line_nocommas.split()
        assert words[0] == "shape"
        shapetype = words[1]
        if shapetype == "sphere":
            self.stru.pdffit["spdiameter"] = float(words[2])
        elif shapetype == "stepcut":
            self.stru.pdffit["stepcut"] = float(words[2])
        else:
            emsg = "Invalid type of particle shape correction %r" % shapetype
            raise StructureFormatError(emsg)
        return


# End of class P_pdffit

# Routines -------------------------------------------------------------------


def getParser():
    """Return new `parser` object for PDFfit format.

    Returns
    -------
    P_pdffit
        Instance of `P_pdffit`.
    """
    return P_pdffit()

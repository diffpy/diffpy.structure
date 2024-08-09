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

"""Basic parser for PDB structure format.

Note
----
References:
    https://www.wwpdb.org/documentation/file-format-content/format23/v2.3.html
    https://www.wwpdb.org/documentation/file-format-content/format30/index.html
"""

import sys

import numpy
from numpy import pi

from diffpy.structure import Structure
from diffpy.structure.parsers import StructureParser
from diffpy.structure.structureerrors import StructureFormatError


class P_pdb(StructureParser):
    """Simple parser for PDB format.

    The parser understands following PDB records: `TITLE, CRYST1, SCALE1,
    SCALE2, SCALE3, ATOM, SIGATM, ANISOU, SIGUIJ, TER, HETATM, END`.

    Attributes
    ----------
    format : str
        Format name, default "pdb".
    """

    # Static data members
    orderOfRecords = [
        "HEADER",
        "OBSLTE",
        "TITLE",
        "CAVEAT",
        "COMPND",
        "SOURCE",
        "KEYWDS",
        "EXPDTA",
        "AUTHOR",
        "REVDAT",
        "SPRSDE",
        "JRNL",
        "REMARK",
        "REMARK",
        "REMARK",
        "REMARK",
        "DBREF",
        "SEQADV",
        "SEQRES",
        "MODRES",
        "HET",
        "HETNAM",
        "HETSYN",
        "FORMUL",
        "HELIX",
        "SHEET",
        "TURN",
        "SSBOND",
        "LINK",
        "HYDBND",
        "SLTBRG",
        "CISPEP",
        "SITE",
        "CRYST1",
        "ORIGX1",
        "ORIGX2",
        "ORIGX3",
        "SCALE1",
        "SCALE2",
        "SCALE3",
        "MTRIX1",
        "MTRIX2",
        "MTRIX3",
        "TVECT",
        "MODEL",
        "ATOM",
        "SIGATM",
        "ANISOU",
        "SIGUIJ",
        "TER",
        "HETATM",
        "ENDMDL",
        "CONECT",
        "MASTER",
        "END",
    ]
    """list: Ordered list of PDB record labels."""

    validRecords = dict.fromkeys(orderOfRecords)
    """dict: Dictionary of PDB record labels."""

    def __init__(self):
        StructureParser.__init__(self)
        self.format = "pdb"
        return

    def parseLines(self, lines):
        """Parse list of lines in PDB format.

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
            Invalid PDB record.
        """
        try:
            stru = Structure()
            scale = numpy.identity(3, dtype=float)
            scaleU = numpy.zeros(3, dtype=float)
            p_nl = 0
            for line in lines:
                p_nl += 1
                # skip blank lines
                if not line.strip():
                    continue
                # make sure line has 80 characters
                if len(line) < 80:
                    line = "%-80s" % line
                words = line.split()
                record = words[0]
                if record == "TITLE":
                    continuation = line[8:10]
                    if continuation.strip():
                        stru.title += line[10:].rstrip()
                    else:
                        stru.title = line[10:].rstrip()
                elif record == "CRYST1":
                    a = float(line[7:15])
                    b = float(line[15:24])
                    c = float(line[24:33])
                    alpha = float(line[33:40])
                    beta = float(line[40:47])
                    gamma = float(line[47:54])
                    stru.lattice.setLatPar(a, b, c, alpha, beta, gamma)
                    scale = numpy.transpose(stru.lattice.recbase)
                elif record == "SCALE1":
                    sc = numpy.zeros((3, 3), dtype=float)
                    sc[0, :] = [float(x) for x in line[10:40].split()]
                    scaleU[0] = float(line[45:55])
                elif record == "SCALE2":
                    sc[1, :] = [float(x) for x in line[10:40].split()]
                    scaleU[1] = float(line[45:55])
                elif record == "SCALE3":
                    sc[2, :] = [float(x) for x in line[10:40].split()]
                    scaleU[2] = float(line[45:55])
                    base = numpy.transpose(numpy.linalg.inv(sc))
                    abcABGcryst = numpy.array(stru.lattice.abcABG())
                    stru.lattice.setLatBase(base)
                    abcABGscale = numpy.array(stru.lattice.abcABG())
                    reldiff = numpy.fabs(1.0 - abcABGscale / abcABGcryst)
                    if not numpy.all(reldiff < 1.0e-4):
                        emsg = "%d: " % p_nl + "SCALE and CRYST1 are not consistent."
                        raise StructureFormatError(emsg)
                    if numpy.any(scaleU != 0.0):
                        emsg = "Origin offset not yet implemented."
                        raise NotImplementedError(emsg)
                elif record in ("ATOM", "HETATM"):
                    name = line[12:16].strip()
                    rc = [float(x) for x in line[30:54].split()]
                    try:
                        occupancy = float(line[54:60])
                    except ValueError:
                        occupancy = 1.0
                    try:
                        B = float(line[60:66])
                        uiso = B / (8 * pi**2)
                    except ValueError:
                        uiso = 0.0
                    element = line[76:78].strip()
                    if element == "":
                        # get element from the first 2 characters of name
                        element = line[12:14].strip()
                        element = element[0].upper() + element[1:].lower()
                    stru.addNewAtom(element, occupancy=occupancy, label=name)
                    last_atom = stru.getLastAtom()
                    last_atom.xyz_cartn = rc
                    last_atom.Uisoequiv = uiso
                elif record == "SIGATM":
                    sigrc = [float(x) for x in line[30:54].split()]
                    sigxyz = numpy.dot(scale, sigrc)
                    try:
                        sigo = float(line[54:60])
                    except ValueError:
                        sigo = 0.0
                    try:
                        sigB = float(line[60:66])
                        sigU = numpy.identity(3) * sigB / (8 * pi**2)
                    except ValueError:
                        sigU = numpy.zeros((3, 3), dtype=float)
                    last_atom.sigxyz = sigxyz
                    last_atom.sigo = sigo
                    last_atom.sigU = sigU
                elif record == "ANISOU":
                    last_atom.anisotropy = True
                    Uij = [float(x) * 1.0e-4 for x in line[28:70].split()]
                    Ua = last_atom.U
                    for i in range(3):
                        Ua[i, i] = Uij[i]
                    Ua[0, 1] = Ua[1, 0] = Uij[3]
                    Ua[0, 2] = Ua[2, 0] = Uij[4]
                    Ua[1, 2] = Ua[2, 1] = Uij[5]
                elif record == "SIGUIJ":
                    sigUij = [float(x) * 1.0e-4 for x in line[28:70].split()]
                    for i in range(3):
                        last_atom.sigU[i, i] = sigUij[i]
                    last_atom.sigU[0, 1] = last_atom.sigU[1, 0] = sigUij[3]
                    last_atom.sigU[0, 2] = last_atom.sigU[2, 0] = sigUij[4]
                    last_atom.sigU[1, 2] = last_atom.sigU[2, 1] = sigUij[5]
                elif record in P_pdb.validRecords:
                    pass
                else:
                    emsg = "%d: invalid record name '%r'" % (p_nl, record)
                    raise StructureFormatError(emsg)
        except (ValueError, IndexError):
            emsg = "%d: invalid PDB record" % p_nl
            exc_type, exc_value, exc_traceback = sys.exc_info()
            e = StructureFormatError(emsg)
            raise e.with_traceback(exc_traceback)
        return stru

    def titleLines(self, stru):
        """Build lines corresponding to `TITLE` record."""
        lines = []
        title = stru.title
        while title != "":
            stop = len(title)
            # maximum length of title record is 60
            if stop > 60:
                stop = title.rfind(" ", 10, 60)
                if stop < 0:
                    stop = 60
            if len(lines) == 0:
                continuation = "  "
            else:
                continuation = "%2i" % (len(lines) + 1)
            lines.append("%-80s" % ("TITLE   " + continuation + title[0:stop]))
            title = title[stop:]
        return lines

    def cryst1Lines(self, stru):
        """Build lines corresponding to `CRYST1` record."""
        lines = []
        latpar = (
            stru.lattice.a,
            stru.lattice.b,
            stru.lattice.c,
            stru.lattice.alpha,
            stru.lattice.beta,
            stru.lattice.gamma,
        )
        if latpar != (1.0, 1.0, 1.0, 90.0, 90.0, 90.0):
            line = "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f" % latpar
            lines.append("%-80s" % line)
        return lines

    def atomLines(self, stru, idx):
        """Build `ATOM` records and possibly `SIGATM`, `ANISOU` or `SIGUIJ` records
        for `structure` stru `atom` number aidx.
        """
        lines = []
        a = stru[idx]
        ad = a.__dict__
        rc = a.xyz_cartn
        B = a.Bisoequiv
        atomline = (
            "ATOM  "  # 1-6
            + "%(serial)5i "  # 7-11, 12
            + "%(name)-4s"  # 13-16
            + "%(altLoc)c"  # 17
            + "%(resName)-3s "  # 18-20, 21
            + "%(chainID)c"  # 22
            + "%(resSeq)4i"  # 23-26
            + "%(iCode)c   "  # 27, 28-30
            + "%(x)8.3f%(y)8.3f%(z)8.3f"  # 31-54
            + "%(occupancy)6.2f"  # 55-60
            + "%(tempFactor)6.2f      "  # 61-66, 67-72
            + "%(segID)-4s"  # 73-76
            + "%(element)2s"  # 77-78
            + "%(charge)-2s"  # 79-80
        ) % {
            "serial": idx + 1,
            "name": a.label or a.element,
            "altLoc": " ",
            "resName": "",
            "chainID": " ",
            "resSeq": 1,
            "iCode": " ",
            "x": rc[0],
            "y": rc[1],
            "z": rc[2],
            "occupancy": a.occupancy,
            "tempFactor": B,
            "segID": "",
            "element": a.element,
            "charge": "",
        }
        lines.append(atomline)
        isotropic = numpy.all(a.U == a.U[0, 0] * numpy.identity(3))
        if not isotropic:
            mid = " %7i%7i%7i%7i%7i%7i  " % tuple(
                numpy.around(1e4 * numpy.array([a.U[0, 0], a.U[1, 1], a.U[2, 2], a.U[0, 1], a.U[0, 2], a.U[1, 2]]))
            )
            line = "ANISOU" + atomline[6:27] + mid + atomline[72:80]
            lines.append(line)
        # default values of standard deviations
        d_sigxyz = numpy.zeros(3, dtype=float)
        d_sigo = 0.0
        d_sigU = numpy.zeros((3, 3), dtype=float)
        sigxyz = ad.get("sigxyz", d_sigxyz)
        sigo = [ad.get("sigo", d_sigo)]
        sigU = ad.get("sigU", d_sigU)
        sigB = [8 * pi**2 * numpy.average([sigU[i, i] for i in range(3)])]
        sigmas = numpy.concatenate((sigxyz, sigo, sigB))
        # no need to print sigmas if they all round to zero
        hassigmas = numpy.any(numpy.fabs(sigmas) >= numpy.array(3 * [5e-4] + 2 * [5e-3])) or numpy.any(
            numpy.fabs(sigU) > 5.0e-5
        )
        if hassigmas:
            mid = "   %8.3f%8.3f%8.3f%6.2f%6.2f      " % tuple(sigmas)
            line = "SIGATM" + atomline[6:27] + mid + atomline[72:80]
            lines.append(line)
            # do we need SIGUIJ record?
            if not numpy.all(sigU == sigU[0, 0] * numpy.identity(3)):
                mid = " %7i%7i%7i%7i%7i%7i  " % tuple(
                    numpy.around(
                        1e4 * numpy.array([sigU[0, 0], sigU[1, 1], sigU[2, 2], sigU[0, 1], sigU[0, 2], sigU[1, 2]])
                    )
                )
                line = "SIGUIJ" + atomline[6:27] + mid + atomline[72:80]
                lines.append(line)
        return lines

    def toLines(self, stru):
        """Convert `Structure` stru to a list of lines in PDB format.

        Parameters
        ----------
        stru : Structure
            Structure to be converted.

        Returns
        -------
        list of str
            List of lines in PDB format.
        """
        lines = []
        lines.extend(self.titleLines(stru))
        lines.extend(self.cryst1Lines(stru))
        for idx in range(len(stru)):
            lines.extend(self.atomLines(stru, idx))
        line = (
            "TER   "  # 1-6
            + "%(serial)5i      "  # 7-11, 12-17
            + "%(resName)-3s "  # 18-20, 21
            + "%(chainID)c"  # 22
            + "%(resSeq)4i"  # 23-26
            + "%(iCode)c"  # 27
            + "%(blank)53s"  # 28-80
        ) % {"serial": len(stru) + 1, "resName": "", "chainID": " ", "resSeq": 1, "iCode": " ", "blank": " "}
        lines.append(line)
        lines.append("%-80s" % "END")
        return lines


# End of class P_pdb

# Routines -------------------------------------------------------------------


def getParser():
    """Return new `parser` object for PDB format.

    Returns
    -------
    P_pdb
        Instance of `P_pdb`.
    """
    return P_pdb()

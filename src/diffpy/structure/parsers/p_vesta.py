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
"""Parser for VESTA format used by VESTA (Visualization for Electronic
and Structural Analysis).

This module replaces the AtomEye XCFG parser (P_xcfg). The XCFG parser and
all its original attributes are preserved for backward compatibility.
VESTA is the actively maintained successor viewer.

Attributes
----------
AtomicMass : dict
    Dictionary of atomic masses for elements.
"""

import sys

import numpy

from diffpy.structure import Structure
from diffpy.structure.parsers import StructureParser
from diffpy.structure.structureerrors import StructureFormatError

# Constants ------------------------------------------------------------------

# Atomic Mass of elements
# This can be later when PeriodicTable package becomes available.

AtomicMass = {
    "H": 1.007947,  # 1 H hydrogen 1.007947
    "He": 4.0026022,  # 2 He helium 4.0026022
    "Li": 6.9412,  # 3 Li lithium 6.9412
    "Be": 9.0121823,  # 4 Be beryllium 9.0121823
    "B": 10.8117,  # 5 B boron 10.8117
    "C": 12.01078,  # 6 C carbon 12.01078
    "N": 14.00672,  # 7 N nitrogen 14.00672
    "O": 15.99943,  # 8 O oxygen 15.99943
    "F": 18.99840325,  # 9 F fluorine 18.99840325
    "Ne": 20.17976,  # 10 Ne neon 20.17976
    "Na": 22.9897702,  # 11 Na sodium 22.9897702
    "Mg": 24.30506,  # 12 Mg magnesium 24.30506
    "Al": 26.9815382,  # 13 Al aluminium 26.9815382
    "Si": 28.08553,  # 14 Si silicon 28.08553
    "P": 30.9737612,  # 15 P phosphorus 30.9737612
    "S": 32.0655,  # 16 S sulfur 32.0655
    "Cl": 35.4532,  # 17 Cl chlorine 35.4532
    "Ar": 39.9481,  # 18 Ar argon 39.9481
    "K": 39.09831,  # 19 K potassium 39.09831
    "Ca": 40.0784,  # 20 Ca calcium 40.0784
    "Sc": 44.9559108,  # 21 Sc scandium 44.9559108
    "Ti": 47.8671,  # 22 Ti titanium 47.8671
    "V": 50.94151,  # 23 V vanadium 50.94151
    "Cr": 51.99616,  # 24 Cr chromium 51.99616
    "Mn": 54.9380499,  # 25 Mn manganese 54.9380499
    "Fe": 55.8452,  # 26 Fe iron 55.8452
    "Co": 58.9332009,  # 27 Co cobalt 58.9332009
    "Ni": 58.69342,  # 28 Ni nickel 58.69342
    "Cu": 63.5463,  # 29 Cu copper 63.5463
    "Zn": 65.4094,  # 30 Zn zinc 65.4094
    "Ga": 69.7231,  # 31 Ga gallium 69.7231
    "Ge": 72.641,  # 32 Ge germanium 72.641
    "As": 74.921602,  # 33 As arsenic 74.921602
    "Se": 78.963,  # 34 Se selenium 78.963
    "Br": 79.9041,  # 35 Br bromine 79.9041
    "Kr": 83.7982,  # 36 Kr krypton 83.7982
    "Rb": 85.46783,  # 37 Rb rubidium 85.46783
    "Sr": 87.621,  # 38 Sr strontium 87.621
    "Y": 88.905852,  # 39 Y yttrium 88.905852
    "Zr": 91.2242,  # 40 Zr zirconium 91.2242
    "Nb": 92.906382,  # 41 Nb niobium 92.906382
    "Mo": 95.942,  # 42 Mo molybdenum 95.942
    "Tc": 98.0,  # 43 Tc technetium 98
    "Ru": 101.072,  # 44 Ru ruthenium 101.072
    "Rh": 102.905502,  # 45 Rh rhodium 102.905502
    "Pd": 106.421,  # 46 Pd palladium 106.421
    "Ag": 107.86822,  # 47 Ag silver 107.86822
    "Cd": 112.4118,  # 48 Cd cadmium 112.4118
    "In": 114.8183,  # 49 In indium 114.8183
    "Sn": 118.7107,  # 50 Sn tin 118.7107
    "Sb": 121.7601,  # 51 Sb antimony 121.7601
    "Te": 127.603,  # 52 Te tellurium 127.603
    "I": 126.904473,  # 53 I iodine 126.904473
    "Xe": 131.2936,  # 54 Xe xenon 131.2936
    "Cs": 132.905452,  # 55 Cs caesium 132.905452
    "Ba": 137.3277,  # 56 Ba barium 137.3277
    "La": 138.90552,  # 57 La lanthanum 138.90552
    "Ce": 140.1161,  # 58 Ce cerium 140.1161
    "Pr": 140.907652,  # 59 Pr praseodymium 140.907652
    "Nd": 144.243,  # 60 Nd neodymium 144.243
    "Pm": 145.0,  # 61 Pm promethium 145
    "Sm": 150.363,  # 62 Sm samarium 150.363
    "Eu": 151.9641,  # 63 Eu europium 151.9641
    "Gd": 157.253,  # 64 Gd gadolinium 157.253
    "Tb": 158.925342,  # 65 Tb terbium 158.925342
    "Dy": 162.5001,  # 66 Dy dysprosium 162.5001
    "Ho": 164.930322,  # 67 Ho holmium 164.930322
    "Er": 167.2593,  # 68 Er erbium 167.2593
    "Tm": 168.934212,  # 69 Tm thulium 168.934212
    "Yb": 173.043,  # 70 Yb ytterbium 173.043
    "Lu": 174.9671,  # 71 Lu lutetium 174.9671
    "Hf": 178.492,  # 72 Hf hafnium 178.492
    "Ta": 180.94791,  # 73 Ta tantalum 180.94791
    "W": 183.841,  # 74 W tungsten 183.841
    "Re": 186.2071,  # 75 Re rhenium 186.2071
    "Os": 190.233,  # 76 Os osmium 190.233
    "Ir": 192.2173,  # 77 Ir iridium 192.2173
    "Pt": 195.0782,  # 78 Pt platinum 195.0782
    "Au": 196.966552,  # 79 Au gold 196.966552
    "Hg": 200.592,  # 80 Hg mercury 200.592
    "Tl": 204.38332,  # 81 Tl thallium 204.38332
    "Pb": 207.21,  # 82 Pb lead 207.21
    "Bi": 208.980382,  # 83 Bi bismuth 208.980382
    "Po": 209.0,  # 84 Po polonium 209
    "At": 210.0,  # 85 At astatine 210
    "Rn": 222.0,  # 86 Rn radon 222
    "Fr": 223.0,  # 87 Fr francium 223
    "Ra": 226.0,  # 88 Ra radium 226
    "Ac": 227.0,  # 89 Ac actinium 227
    "Th": 232.03811,  # 90 Th thorium 232.03811
    "Pa": 231.035882,  # 91 Pa protactinium 231.035882
    "U": 238.028913,  # 92 U uranium 238.028913
    "Np": 237.0,  # 93 Np neptunium 237
    "Pu": 244.0,  # 94 Pu plutonium 244
    "Am": 243.0,  # 95 Am americium 243
    "Cm": 247.0,  # 96 Cm curium 247
    "Bk": 247.0,  # 97 Bk berkelium 247
    "Cf": 251.0,  # 98 Cf californium 251
    "Es": 252.0,  # 99 Es einsteinium 252
    "Fm": 257.0,  # 100 Fm fermium 257
    "Md": 258.0,  # 101 Md mendelevium 258
    "No": 259.0,  # 102 No nobelium 259
    "Lr": 262.0,  # 103 Lr lawrencium 262
    "Rf": 261.0,  # 104 Rf rutherfordium 261
    "Db": 262.0,  # 105 Db dubnium 262
    "Sg": 266.0,  # 106 Sg seaborgium 266
    "Bh": 264.0,  # 107 Bh bohrium 264
    "Hs": 277.0,  # 108 Hs hassium 277
    "Mt": 268.0,  # 109 Mt meitnerium 268
    "Ds": 281.0,  # 110 Ds darmstadtium 281
    "Rg": 272.0,  # 111 Rg roentgenium 272
}


class P_vesta(StructureParser):
    """Parser for VESTA native structure format (.vesta).

    VESTA (Visualization for Electronic and Structural Analysis) is the
    actively maintained successor to AtomEye. This parser writes the
    native VESTA format understood by VESTA 3.x and later.

    Attributes
    ----------
    format : str
        Format name, default "vesta".

    Notes
    -----
    The ``cluster_boundary`` attribute is retained from the original
    AtomEye/XCFG parser for API compatibility; it is not used by VESTA
    because VESTA handles periodicity natively.
    """

    cluster_boundary = 2
    """int: Width of boundary around corners of non-periodic cluster.
    Retained from the original AtomEye/XCFG parser for API compatibility.
    VESTA handles periodicity natively so this value has no effect on output.
    """

    def __init__(self):
        StructureParser.__init__(self)
        self.format = "vesta"
        return

    def parse_lines(self, lines):
        """Parse list of lines in VESTA format.

        Reads the ``STRUC``, ``ATOMT``, and ``COORD`` sections of a
        ``.vesta`` file to reconstruct a :class:`~diffpy.structure.Structure`.

        Parameters
        ----------
        lines : list of str
            Lines of a VESTA format file.

        Returns
        -------
        Structure
            Parsed structure instance.

        Raises
        ------
        StructureFormatError
            When the file does not conform to the VESTA format.
        """
        stru = Structure()
        p_nl = 0

        # Strip trailing blank lines for a clean iteration boundary.
        stop = len(lines)
        for line in reversed(lines):
            if line.strip():
                break
            stop -= 1
        ilines = iter(lines[:stop])

        try:
            # Lattice parameters parsed from STRUC block:
            #   a b c alpha beta gamma
            latt_abc = None
            latt_abg = None
            atom_types = {}

            # Raw fractional coordinates collected from COORD block:
            #   list of (atom_type_index, x, y, z, occupancy)
            raw_coords = []

            section = None  # tracks current block keyword

            for line in ilines:
                p_nl += 1
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue

                # Detect section transitions.
                upper = stripped.split()[0].upper()
                if upper in (
                    "CRYSTAL",
                    "TITLE",
                    "GROUP",
                    "STRUC",
                    "ATOMT",
                    "COORD",
                    "BOUND",
                    "SBOND",
                    "VECTR",
                    "VECTS",
                    "STYLE",
                    "SCENE",
                    "EOF",
                ):
                    section = upper
                    continue

                # ---- STRUC section: lattice parameters -----------------
                if section == "STRUC":
                    words = stripped.split()
                    # First data line: a b c  alpha beta gamma  space_group
                    if latt_abc is None and len(words) >= 6:
                        try:
                            latt_abc = [float(w) for w in words[:3]]
                            latt_abg = [float(w) for w in words[3:6]]
                        except ValueError:
                            pass
                    continue

                # ---- ATOMT section: atom-type definitions ---------------
                if section == "ATOMT":
                    # Format: index  Symbol  radius  r  g  b  ...
                    words = stripped.split()
                    if len(words) >= 2:
                        try:
                            idx = int(words[0])
                            symbol = words[1]
                            atom_types[idx] = symbol
                        except ValueError:
                            pass
                    continue

                # ---- COORD section: atomic coordinates -----------------
                if section == "COORD":
                    # Format: seq  type_index  x  y  z  occupancy  ...
                    words = stripped.split()
                    if len(words) >= 6:
                        try:
                            type_idx = int(words[1])
                            x, y, z = float(words[2]), float(words[3]), float(words[4])
                            occ = float(words[5])
                            raw_coords.append((type_idx, x, y, z, occ))
                        except ValueError:
                            pass
                    continue
            if latt_abc is None:
                emsg = "VESTA file is missing STRUC lattice parameters"
                raise StructureFormatError(emsg)

            stru.lattice.setLatPar(
                a=latt_abc[0],
                b=latt_abc[1],
                c=latt_abc[2],
                alpha=latt_abg[0],
                beta=latt_abg[1],
                gamma=latt_abg[2],
            )

            for type_idx, x, y, z, occ in raw_coords:
                element = atom_types.get(type_idx, "X")
                stru.add_new_atom(element, xyz=[x, y, z])
                stru[-1].occupancy = occ

        except (ValueError, IndexError):
            emsg = "%d: file is not in VESTA format" % p_nl
            exc_type, exc_value, exc_traceback = sys.exc_info()
            e = StructureFormatError(emsg)
            raise e.with_traceback(exc_traceback)

        return stru

    def to_lines(self, stru):
        """Convert Structure *stru* to a list of lines in VESTA format.

        Produces a ``.vesta`` file readable by VESTA 3.x and later,
        containing ``STRUC``, ``ATOMT``, and ``COORD`` sections derived
        from the structure's lattice and atomic data.

        Parameters
        ----------
        stru : Structure
            Structure to be converted.

        Returns
        -------
        list of str
            Lines of a VESTA format file.

        Raises
        ------
        StructureFormatError
            Cannot convert empty structure to VESTA format.
        """
        if len(stru) == 0:
            emsg = "cannot convert empty structure to VESTA format"
            raise StructureFormatError(emsg)

        lines = []
        lines.append("#VESTA_FORMAT_VERSION 3.5.0")
        lines.append("")
        lines.append("CRYSTAL")
        lines.append("")
        lines.append("TITLE")
        title = getattr(stru, "title", "") or "Structure"
        lines.append(title)
        lines.append("")
        latt = stru.lattice
        a, b, c, alpha, beta, gamma = latt.cell_parms()
        lines.append("STRUC")
        # Line 1: a  b  c  alpha  beta  gamma  space_group_number
        lines.append("  %.8g  %.8g  %.8g  %.8g  %.8g  %.8g  1" % (a, b, c, alpha, beta, gamma))
        # Line 2: origin shift (0 0 0) followed by space-group symbol placeholder
        lines.append("  0.000000  0.000000  0.000000")
        lines.append("")
        element_order = []
        seen = set()
        for a_obj in stru:
            el = a_obj.element
            if el not in seen:
                seen.add(el)
                element_order.append(el)
        type_index = {el: i + 1 for i, el in enumerate(element_order)}
        lines.append("ATOMT")
        for el in element_order:
            idx = type_index[el]
            # Default ball radius 0.5; placeholder RGB 1.0 1.0 1.0.
            lines.append("  %d  %s  %.4f  1.0000  1.0000  1.0000  204" % (idx, el, 0.5))
        lines.append("")
        lines.append("COORD")
        for seq, a_obj in enumerate(stru, start=1):
            el = a_obj.element
            tidx = type_index[el]
            x, y, z = a_obj.xyz
            occ = getattr(a_obj, "occupancy", 1.0)
            # Isotropic displacement parameter (Uiso), defaulting to 0.
            uiso = _get_uiso(a_obj)
            lines.append("  %d  %d  %.8g  %.8g  %.8g  %.4f  %.4f" % (seq, tidx, x, y, z, occ, uiso))
        lines.append("  0 0 0 0 0")
        lines.append("")
        lines.append("BOUND")
        lines.append("  0.0  1.0  0.0  1.0  0.0  1.0")
        lines.append("  0  0  0  0  0")
        lines.append("")
        lines.append("EOF")
        return lines
# End of class P_vesta

from diffpy.structure.parsers.P_xcfg import P_xcfg  # noqa: E402, F401

# Routines -------------------------------------------------------------------


def get_parser():
    """Return new parser object for VESTA format.

    Returns
    -------
    P_vesta
        Instance of :class:`P_vesta`.
    """
    return P_vesta()


# Local Helpers --------------------------------------------------------------


def _get_uiso(a):
    """Return isotropic displacement parameter for atom *a*.

    Tries ``Uisoequiv`` first, then falls back to the mean of the
    diagonal of the anisotropic U tensor, then to zero.

    Parameters
    ----------
    a : Atom
        Atom instance.

    Returns
    -------
    float
        Isotropic U value in Å².
    """
    if hasattr(a, "Uisoequiv"):
        return float(a.Uisoequiv)
    try:
        return float(numpy.trace(a.U) / 3.0)
    except Exception:
        return 0.0


def _assign_auxiliaries(a, fields, auxiliaries, no_velocity):
    """Assign auxiliary properties for an
    :class:`~diffpy.structure.Atom` object.

    Retained from the original AtomEye/XCFG parser for backward
    compatibility with code that calls this helper directly.

    Parameters
    ----------
    a : Atom
        The Atom instance for which auxiliary properties need to be set.
    fields : list
        Floating-point values for the current row of the processed file.
    auxiliaries : dict
        Dictionary of zero-based indices and names of auxiliary properties.
    no_velocity : bool
        When ``False``, set atom velocity ``a.v`` to ``fields[3:6]``.
        Use ``fields[3:6]`` for auxiliary values otherwise.
    """
    if not no_velocity:
        a.v = numpy.asarray(fields[3:6], dtype=float)
    auxfirst = 3 if no_velocity else 6
    for i, prop in auxiliaries.items():
        value = fields[auxfirst + i]
        if prop == "Uiso":
            a.Uisoequiv = value
        elif prop == "Biso":
            a.Bisoequiv = value
        elif prop[0] in "BU" and all(d in "123" for d in prop[1:]):
            nm = prop if prop[1] <= prop[2] else prop[0] + prop[2] + prop[1]
            a.anisotropy = True
            setattr(a, nm, value)
        else:
            setattr(a, prop, value)
    return

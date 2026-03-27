#!/usr/bin/env python
##############################################################################
#
# diffpy.structure  by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2026 University of California, Santa Barbara.
#                   All rights reserved.
#
# File coded by:    Simon J. L. Billinge, Rundong Hua
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

import re
import sys

import numpy

from diffpy.structure import Structure
from diffpy.structure.parsers import StructureParser
from diffpy.structure.parsers.p_xcfg import AtomicMass
from diffpy.structure.structureerrors import StructureFormatError


# Constants ------------------------------------------------------------------
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
                    # Format: index  Symbol  radius  r  g  b  style  # mass=<value>
                    words = stripped.split()
                    if len(words) >= 2:
                        try:
                            idx = int(words[0])
                            symbol = words[1]
                            atom_types[idx] = {"symbol": symbol, "mass": None}
                            # Recover mass from the trailing comment if present.
                            mass_match = re.search(r"#\s*mass\s*=\s*([0-9.eE+\-]+)", stripped)
                            if mass_match:
                                atom_types[idx]["mass"] = float(mass_match.group(1))
                            else:
                                # Fall back to the built-in lookup table.
                                atom_types[idx]["mass"] = AtomicMass.get(symbol, 0.0)
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
                type_info = atom_types.get(type_idx, {"symbol": "X", "mass": 0.0})
                element = type_info["symbol"]
                mass = type_info["mass"]
                stru.add_new_atom(element, xyz=[x, y, z])
                stru[-1].occupancy = occ
                if mass is None:
                    mass = AtomicMass.get(element, 0.0)
                stru[-1].mass = mass

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
            mass = AtomicMass.get(el, 0.0)
            lines.append("  %d  %s  %.4f  1.0000  1.0000  1.0000  204  # mass=%.7g" % (idx, el, 0.5, mass))
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

#!/usr/bin/env python
##############################################################################
#
# diffpy.Structure  by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2007 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    C. L. Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################

"""Parser for Bruce Ravel's Atoms structure format
"""

from diffpy.Structure import Lattice, Atom, StructureFormatError
from diffpy.Structure.bratomsstructure import BRAtomsStructure
from diffpy.Structure.Parsers import StructureParser

class P_bratoms(StructureParser):
    """Parser for Bruce Ravel's Atoms structure format.
    """

    plist = ["a", "b", "c", "alpha", "beta", "gamma"]

    def __init__(self):
        StructureParser.__init__(self)
        self.format = "bratoms"
        return

    def parseLines(self, lines):
        """Parse list of lines in atoms format.

        Return Structure object or raise StructureFormatError.
        """

        comlist = ["#", "%", "!", "*"]
        atoms = []
        title = ""
        anext = False
        structure = BRAtomsStructure()
        meta = structure.bratoms
        pdict = dict.fromkeys(self.plist)


        # Count the lines
        ln = 0
        try:

            for line in lines:
                ln += 1

                # Strip comments from the line
                for c in comlist:
                    idx = line.find(c)
                    if idx != -1:
                        line = line[:idx]

                # Move on if there is not a line
                if not line: continue

                # Move on if there was only white space in the line
                sline = line.split()
                if not sline: continue

                # Check if we have atoms following
                if sline[0].startswith("atom"):
                    anext = True
                    continue

                # Check for title
                if sline[0].startswith("title"):
                    if title: title += "\n"
                    title += line[5:]
                    continue

                # Get rid of pesky "=" and "," signs
                while "=" in sline: sline.remove("=")
                while "," in sline: sline.remove(",")

                # space group
                if sline and sline[0].startswith("space"):
                    meta["space"] = line[5:].strip()
                    continue

                # output
                if sline and sline[0].startswith("output"):
                    meta["output"] = line[6:].strip()
                    continue

                # shift
                if sline and sline[0].startswith("shift"):
                    meta["shift"] = line[5:].strip()
                    continue

                # Check for other metadata
                while sline and sline[0].strip() in meta:
                    key = sline.pop(0).strip()
                    if key == "central": key = "core"
                    meta[key] = sline.pop(0).strip()

                # Check for lattice information.
                while sline and sline[0].strip() in self.plist:
                    key = sline.pop(0).strip()
                    pdict[key] = float(sline.pop(0))

                # Check for atom information
                if sline and anext:

                    elraw = sline.pop(0).strip()
                    el = elraw[:1].upper() + elraw[1:].lower()
                    x = float(sline.pop(0))
                    y = float(sline.pop(0))
                    z = float(sline.pop(0))

                    tag = ""
                    if sline:
                        tag = sline.pop(0).strip()
                    occ = 1.0
                    if sline:
                        occ = float(sline.pop(0))

                    a = Atom(atype = el,
                        xyz = [x, y, z],
                        label = tag,
                        occupancy = occ)

                    atoms.append(a)

        except (ValueError, IndexError):
            emsg = "%d: file is not in Atoms format" % ln
            raise StructureFormatError(emsg)

        # Make sure we have atoms.
        if len(atoms) == 0:
            raise StructureFormatError("File contains no atoms")

        # Make sure we have unit cell parameters
        if pdict["a"] is None:
            emsg = "Missing definition of cell parameter"
            raise StructureFormatError(emsg)

        # Fill in optional information if it was missing.
        if pdict["alpha"] is None:
            pdict["alpha"] = 90.0
        if pdict["beta"] is None:
            pdict["beta"] = pdict["alpha"]
        if pdict["gamma"] is None:
            pdict["gamma"] = pdict["alpha"]
        if pdict["b"] is None:
            pdict["b"] = pdict["a"]
        if pdict["c"] is None:
            pdict["c"] = pdict["a"]

        if meta['core'] is None:
            meta['core'] = atoms[0].element

        lat = Lattice(**pdict)
        structure.title = title
        structure.lattice = lat
        structure.extend(atoms)
        return structure
    # End of parseLines

    def toLines(self, stru):
        """Convert Structure stru to a list of lines in PDFfit format.

        Return list of strings.
        """
        lat = stru.lattice

        # Check the metadata for consistency
        meta = getattr(stru, 'bratoms', {})

        # Make the lines
        lines = []

        # title
        titles = stru.title.split("\n")
        for t in titles:
            lines.append("title = %s" % t)

        # a, b, c
        lines.append("a = %f    b = %f    c = %f" % (lat.a, lat.b, lat.c))
        # alpha, beta, gamma
        lines.append("alpha = %f    beta = %f    gamma = %f" % \
                (lat.alpha, lat.beta, lat.gamma))

        # space
        lines.append("space = %s" % meta.get("space", "P 1"))

        # edge
        lines.append("edge = %s" % meta.get("edge", "K"))

        # core
        tag = meta.get("core") or stru[0].label or stru[0].element.title()
        lines.append("core = %s" % tag)

        # rmax
        lines.append("rmax = %s" % meta.get("rmax", 6.0))

        # Other keys
        if meta.get('shift'):
            lines.append("shift = %s" % meta["shift"])
        if meta.get('output'):
            lines.append("output = %s" % meta["output"])
        if meta.get('nitrogen'):
            lines.append("nitrogen = %s" % meta["nitrogen"])
        if meta.get('argon'):
            lines.append("argon = %s" % meta["argon"])
        if meta.get('krypton'):
            lines.append("krypton = %s" % meta["krypton"])

        # Add atoms
        lines.append("atoms")
        lines.append("%-8s%-10s%-10s%-10s%-10s%s" %
                ("!", "x", "y", "z", "tag", "occ"))
        for a in stru:
            el = a.element.title()
            name = a.label or el
            x, y, z = a.xyz
            occ = a.occupancy
            lines.append("%-7s %-9f %-9f %-9f %-9s %f" % \
                    (el, x, y, z, name, occ))

        return lines
    # End of toLines

# End of class P_pdffit

# Routines

def getParser():
    return P_bratoms()

# End of file

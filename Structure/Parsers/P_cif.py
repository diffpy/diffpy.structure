########################################################################
#
# Structure         by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2006 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
########################################################################

"""Parser for basic CIF file format

http://www.iucr.org/iucr-top/cif/home.html
"""

__id__ = "$Id$"

import sys
import time
import numpy
import numpy.linalg as numalg
from numpy import pi

from Structure.Structure import Structure
from Structure.Lattice import Lattice
from Structure.Atom import Atom
from StructureParser import StructureParser
from Structure.exceptions import InvalidStructureFormat

class Parser(StructureParser):
    """Parser --> StructureParser subclass for basic parser for CIF format"""

    def __init__(self):
        self.format = "cif"
        return

#   parse Lines is not implemented, we probably want to use babel here
#   def parseLines(self, lines):
#       """Parse list of lines in CIF format.
#
#       Return Structure instance or raise InvalidStructureFormat.
#       """
#   # End of parseLines

    def toLines(self, stru):
        """Convert Structure stru to a list of lines in basic CIF format.

        Return list of strings.
        """
        lines = []
        # may be replaced with filtered Structure.title
        # for now, we can add the title as a comment
        if stru.title.strip() != "":
            title_lines = stru.title.split('\n')
            lines.extend([ "# " + line.strip() for line in title_lines ])
            lines.append("")
        lines.append("data_3D")
        iso_date =  "%04i-%02i-%02i" % time.gmtime()[:3]
        lines.extend([
            "%-31s %s" % ("_audit_creation_date", iso_date),
            "%-31s %s" % ("_audit_creation_method", "P_cif.py"),
            "",
            "%-31s %s" % ("_symmetry_space_group_name_H-M", "'P1'"),
            "%-31s %s" % ("_symmetry_Int_Tables_number", "1"),
            "%-31s %s" % ("_symmetry_cell_setting", "triclinic"),
            "" ])
        # there should be no need to specify equivalent positions for P1
        # _symmetry_equiv_posi_as_xyz x,y,z
        lines.extend([
            "%-31s %.6g" % ("_cell_length_a", stru.lattice.a),
            "%-31s %.6g" % ("_cell_length_b", stru.lattice.b),
            "%-31s %.6g" % ("_cell_length_c", stru.lattice.c),
            "%-31s %.6g" % ("_cell_angle_alpha", stru.lattice.alpha),
            "%-31s %.6g" % ("_cell_angle_beta", stru.lattice.beta),
            "%-31s %.6g" % ("_cell_angle_gamma", stru.lattice.gamma),
            "" ])
        # build a list of site labels and adp (displacement factor) types
        element_count = {}
        a_site_label = []
        a_adp_type = []
        for a in stru:
            cnt = element_count[a.element] = element_count.get(a.element,0)+1
            a_site_label.append( "%s%i" % (a.element, cnt) )
            if numpy.all(a.U == a.U[0,0]*numpy.identity(3)):
                a_adp_type.append("Uiso")
            else:
                a_adp_type.append("Uani")
        # list all atoms
        lines.extend([
            "loop_",
            "  _atom_site_label",
            "  _atom_site_type_symbol",
            "  _atom_site_fract_x",
            "  _atom_site_fract_y",
            "  _atom_site_fract_z",
            "  _atom_site_U_iso_or_equiv",
            "  _atom_site_adp_type",
            "  _atom_site_occupancy" ])
        for i in range(len(stru)):
            a = stru[i]
            line = "  %-5s %-3s %11.6f %11.6f %11.6f %11.6f %-5s %.4f" % (
                    a_site_label[i], a.element, a.xyz[0], a.xyz[1], a.xyz[2],
                    a.Uiso(), a_adp_type[i], a.occupancy  )
            lines.append(line)
        # find anisotropic atoms
        idx_aniso = [ i for i in range(len(stru)) if a_adp_type[i] != "Uiso" ]
        if idx_aniso != []:
            lines.extend([
                "loop_",
                "  _atom_site_aniso_label",
                "  _atom_site_aniso_U_11",
                "  _atom_site_aniso_U_22",
                "  _atom_site_aniso_U_33",
                "  _atom_site_aniso_U_12",
                "  _atom_site_aniso_U_13",
                "  _atom_site_aniso_U_23" ])
            for i in idx_aniso:
                a = stru[i]
                line = "  %-5s %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f" % (
                        a_site_label[i], a.U[0,0], a.U[1,1], a.U[2,2],
                        a.U[0,1], a.U[0,2], a.U[1,2] )
                lines.append(line)
        return lines
    # End of toLines

# End of Parser

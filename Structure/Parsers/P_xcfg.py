########################################################################
#
# <PackageName>     by DANSE Diffraction group
#                   Simon J.L. Billinge
#                   Michigan State University
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See COPYRIGHT.txt for copying and usage conditions.
# See LICENSE.txt for license information.
#
########################################################################

"""Parser for extended CFG format used by atomeye"""

__id__ = "$Id$"

import sys
import re
import numpy as num
from Structure.Structure import Structure
from Structure.Lattice import Lattice
from Structure.Atom import Atom
from StructureParser import StructureParser
from Structure.utils import isfloat
from Structure.exceptions import InvalidStructureFormat

class Parser(StructureParser):
    """Parser --> StructureParser subclass for extended CFG format"""

    def __init__(self):
        self.format = "xcfg"
        return

    def parseLines(self, lines):
        """Parse list of lines in PDB format.

        Return Structure object or raise InvalidStructureFormat.
        """
        xcfg_Number_of_particles = None
        xcfg_A = None
        xcfg_H0 = num.zeros((3,3), dtype=float)
        xcfg_H0_set = num.zeros((3,3), dtype=bool)
        xcfg_NO_VELOCITY = False
        xcfg_entry_count = None
        xcfg_auxiliary = []
        p_nl = 0
        p_auxiliary_re = re.compile(r"^auxiliary\[(\d+)\] =")
        p_auxiliary = {}
        try:
            stru = Structure()
            # ignore trailing blank lines
            stop = len(lines)
            while stop>0 and lines[stop-1].strip() == "":
                stop -= 1
            ilines = iter(lines[:stop])
            # read XCFG header
            for line in ilines:
                p_nl += 1
                stripped_line = line.strip()
                # blank lines and lines starting with # are ignored
                if stripped_line == "" or line[0] == '#':
                    continue
                elif xcfg_Number_of_particles is None:
                    if line.find("Number of particles =") != 0:
                        raise InvalidStructureFormat, ("%d: first line must "+
                                "contain 'Number of particles ='") % p_nl
                    xcfg_Number_of_particles = int(line[21:].split(None, 1)[0])
                    p_natoms = xcfg_Number_of_particles
                elif line.find("A =") == 0:
                    xcfg_A = float(line[3:].split(None, 1)[0])
                elif line.find("H0(") == 0:
                    i, j = ( int(line[3])-1 ,  int(line[5])-1 )
                    xcfg_H0[i,j] = float(line[10:].split(None, 1)[0])
                    xcfg_H0_set[i,j] = True
                elif line.find(".NO_VELOCITY.") == 0:
                    xcfg_NO_VELOCITY = True
                elif line.find("entry_count =") == 0:
                    xcfg_entry_count = int(line[13:].split(None, 1)[0])
                elif p_auxiliary_re.match(line):
                    m = p_auxiliary_re.match(line)
                    idx = int(m.group(1))
                    p_auxiliary[idx] = line[m.end():].split(None, 1)[0]
                else:
                    break
            # check header for consistency
            if num.any(xcfg_H0_set == False):
                raise InvalidStructureFormat, \
                        "H0 tensor is not properly defined"
            p_auxnum = len(p_auxiliary) and max(p_auxiliary.keys())+1
            for i in range(p_auxnum):
                if not i in p_auxiliary:
                    p_auxiliary[i] = "aux%d" % i
            sorted_aux_keys = p_auxiliary.keys()
            sorted_aux_keys.sort()
            if p_auxnum != 0:
                stru.xcfg = {
                    'auxiliaries' : [ p_auxiliary[k]
                                      for k in sorted_aux_keys ]
                }
            if 6-3*xcfg_NO_VELOCITY+len(p_auxiliary) != xcfg_entry_count:
                raise InvalidStructureFormat, ("%d: auxiliary fields " +
                        "not consistent with entry_count") % p_nl
            # define proper lattice
            stru.lattice.setLatBase(xcfg_H0)
            # build p_assign_atom function to assign entries to proper fields
            p_exprs = [ "a.xyz[0]=fields[0]",
                        "a.xyz[1]=fields[1]",
                        "a.xyz[2]=fields[2]" ]
            if not xcfg_NO_VELOCITY:
                p_exprs += [  "a.v=num.zeros(3, dtype=float)",
                              "a.v[0]=fields[3]",
                              "a.v[1]=fields[4]",
                              "a.v[2]=fields[5]" ]
            for idx in sorted_aux_keys:
                prop = p_auxiliary[idx]
                col = idx + 6 - 3*xcfg_NO_VELOCITY
                if prop == "Uiso":
                    p_exprs.append("a.U[0,0]=a.U[1,1]=a.U[2,2]=" +
                        "fields[%d]" % col)
                elif re.match(r"^U\d\d$", prop) \
                and 1<=int(prop[1])<=3 and 1<=int(prop[2])<=3 :
                    i, j = int(prop[1])-1, int(prop[2])-1
                    if i==j:
                        p_exprs.append("a.U[%i,%i]=fields[%d]" % (i, j, col) )
                    else:
                        p_exprs.append("a.U[%i,%i]=a.U[%i,%i]=fields[%d]" % \
                                (i, j, j, i, col) )
                else:
                    p_exprs.append( "a.__dict__[%r]=fields[%d]" % \
                            (prop, col) )
            p_assign_expr = "pass; " + "; ".join(p_exprs[3:])
            exec "def p_assign_atom(a, fields) : %s" % p_assign_expr
            # here we are inside data
            p_element = None
            p_nl -= 1
            for line in lines[p_nl:stop]:
                p_nl += 1
                words = line.split()
                if len(words) == 1:
                    # ignore atom mass
                    if isfloat(words[0]):
                        continue
                    else:
                        w = words[0]
                        p_element = w[0].upper() + w[1:].lower()
                elif len(words) == xcfg_entry_count and p_element is not None:
                    fields = [ float(w) for w in words ]
                    a = Atom(p_element, fields[:3])
                    a.xyz *= xcfg_A
                    p_assign_atom(a, fields)
                    stru.append(a)
                else:
                    raise InvalidStructureFormat, "%d: invalid record" % p_nl
            if len(stru) != p_natoms:
                raise InvalidStructureFormat, \
                        "expected %d atoms, read %d" % (p_natoms, len(stru))
        except (ValueError, IndexError):
            exc_type, exc_value, exc_traceback = sys.exc_info()
            raise InvalidStructureFormat, \
                    "%d: file is not in XCFG format" % p_nl, exc_traceback
        return stru
    # End of parseLines

    def toLines(self, stru):
        """Convert Structure stru to a list of lines in XCFG atomeye format.

        Return list of strings.
        """
        from Structure.PeriodicTable import AtomicMass
        if len(stru) == 0:
            raise InvalidStructureFormat, \
                    "cannot convert empty structure to XCFG format"
        lines = []
        lines.append( "Number of particles = %i" % len(stru) )
        # figure out length unit A
        allxyz = num.array([a.xyz for a in stru])
        lo_xyz = num.array([ allxyz[:,i].min() for i in range(3) ])
        hi_xyz = num.array([ allxyz[:,i].max() for i in range(3) ])
        max_range_xyz = (hi_xyz-lo_xyz).max()
        # range of CFG coordinates must be less than 1
        p_A = num.ceil(max_range_xyz + 1.0e-13)
        # atomeye draws rubbish when boxsize is less than 3.5
        hi_ucvect = max([num.sqrt(num.dot(v,v)) for v in stru.lattice.base])
        if hi_ucvect*p_A < 3.5:
            p_A = num.ceil(3.5 / hi_ucvect)
        lines.append( "A = %.8g Angstrom" % p_A )
        # how much do we need to shift the coordinates?
        p_dxyz = num.zeros(3, dtype=float)
        for i in range(3):
            if lo_xyz[i]/p_A < 0.0 or hi_xyz[i]/p_A >= 1.0 \
            or (lo_xyz[i] == hi_xyz[i] and lo_xyz[i] == 0.0) :
                p_dxyz[i] = 0.5 - (hi_xyz[i]+lo_xyz[i])/2.0/p_A
        # H0 tensor
        for i in range(3):
            for j in range(3):
                lines.append( "H0(%i,%i) = %.8g A" % \
                        (i+1, j+1, stru.lattice.base[i,j]) )
        # get out for empty structure
        if len(stru) == 0: return lines
        a_first = stru[0]
        p_NO_VELOCITY = "v" not in a_first.__dict__
        if p_NO_VELOCITY:
            lines.append(".NO_VELOCITY.")
        # build a p_auxiliaries list of (aux_name,atom_expression) tuples
        # if stru came from xcfg file, it would store original auxiliaries in
        # xcfg dictionary
        try:
            p_auxiliaries = [ (aux, "a."+aux)
                for aux in stru.xcfg['auxiliaries'] ]
        except AttributeError:
            p_auxiliaries = []
        # add occupancy if any atom has nonunit occupancy
        for a in stru:
            if a.occupancy != 1.0:
                p_auxiliaries.append( ('occupancy', 'a.occupancy') )
                break
        # add temperature factor with as many terms as needed
        # check whether all temperature factors are zero or isotropic
        p_allUzero = True
        p_allUiso = True
        for a in stru:
            if p_allUzero and num.any(a.U != 0.0):
                p_allUzero = False
            if not num.all(a.U == a.U[0,0]*num.identity(3)):
                p_allUiso = False
                # here p_allUzero must be false
                break
        if p_allUzero:
            pass
        elif p_allUiso:
            p_auxiliaries.append( ('Uiso', 'a.U[0,0]') )
        else:
            p_auxiliaries.extend([ ('U11', 'a.U[0,0]'),
                                   ('U22', 'a.U[1,1]'),
                                   ('U33', 'a.U[2,2]') ])
            # check if there are off-diagonal elements
            allU = num.array([a.U for a in stru])
            if num.any(allU[:,0,1] != 0.0):
                p_auxiliaries.append( ('U12', 'a.U[0,1]') )
            if num.any(allU[:,0,2] != 0.0):
                p_auxiliaries.append( ('U13', 'a.U[0,2]') )
            if num.any(allU[:,1,2] != 0.0):
                p_auxiliaries.append( ('U23', 'a.U[1,2]') )
        # count entries
        p_entry_count = 6 - 3*p_NO_VELOCITY + len(p_auxiliaries)
        lines.append("entry_count = %d" % p_entry_count)
        # add auxiliaries
        for i in range(len(p_auxiliaries)):
            lines.append("auxiliary[%d] = %s [au]" % (i, p_auxiliaries[i][0]))
        # now define p_entry_line function for representing atom properties
        p_exprs = [ "def p_entry_line(a, p_A, p_dxyz):",
                    "    fields = list( a.xyz/p_A+p_dxyz )" ]
        if not p_NO_VELOCITY:
            p_exprs.append( \
                    "    fields += [ a.v[0], a.v[1], a.v[2] ]" )
        p_exprs += ["    fields += [ " +
                         ",".join([e for p,e in p_auxiliaries]) + " ]",
                    "    line = ' '.join([ '%.8g' % x for x in fields ])",
                    "    return line"  ]
        exec "\n".join(p_exprs)
        # we are ready to output atoms:
        lines.append("")
        p_element = None
        for a in stru:
            if a.element != p_element:
                p_element = a.element
                lines.append("%.4f" % AtomicMass.get(p_element, 0.0))
                lines.append(p_element)
            lines.append(p_entry_line(a, p_A, p_dxyz))
        return lines
    # End of toLines

# End of Parser

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

"""Parser for extended CFG format used by atomeye"""

import sys
import re
import numpy

from diffpy.Structure import Structure
from diffpy.Structure import StructureFormatError
from diffpy.Structure.utils import isfloat
from diffpy.Structure.Parsers import StructureParser


##############################################################################
# Constants

# Atomic Mass of elements
# This can be later when PeriodicTable package becomes available.

AtomicMass = {
    "H"  : 1.007947,    # 1 H hydrogen 1.007947
    "He" : 4.0026022,   # 2 He helium 4.0026022
    "Li" : 6.9412,      # 3 Li lithium 6.9412
    "Be" : 9.0121823,   # 4 Be beryllium 9.0121823
    "B"  : 10.8117,     # 5 B boron 10.8117
    "C"  : 12.01078,    # 6 C carbon 12.01078
    "N"  : 14.00672,    # 7 N nitrogen 14.00672
    "O"  : 15.99943,    # 8 O oxygen 15.99943
    "F"  : 18.99840325, # 9 F fluorine 18.99840325
    "Ne" : 20.17976,    # 10 Ne neon 20.17976
    "Na" : 22.9897702,  # 11 Na sodium 22.9897702
    "Mg" : 24.30506,    # 12 Mg magnesium 24.30506
    "Al" : 26.9815382,  # 13 Al aluminium 26.9815382
    "Si" : 28.08553,    # 14 Si silicon 28.08553
    "P"  : 30.9737612,  # 15 P phosphorus 30.9737612
    "S"  : 32.0655,     # 16 S sulfur 32.0655
    "Cl" : 35.4532,     # 17 Cl chlorine 35.4532
    "Ar" : 39.9481,     # 18 Ar argon 39.9481
    "K"  : 39.09831,    # 19 K potassium 39.09831
    "Ca" : 40.0784,     # 20 Ca calcium 40.0784
    "Sc" : 44.9559108,  # 21 Sc scandium 44.9559108
    "Ti" : 47.8671,     # 22 Ti titanium 47.8671
    "V"  : 50.94151,    # 23 V vanadium 50.94151
    "Cr" : 51.99616,    # 24 Cr chromium 51.99616
    "Mn" : 54.9380499,  # 25 Mn manganese 54.9380499
    "Fe" : 55.8452,     # 26 Fe iron 55.8452
    "Co" : 58.9332009,  # 27 Co cobalt 58.9332009
    "Ni" : 58.69342,    # 28 Ni nickel 58.69342
    "Cu" : 63.5463,     # 29 Cu copper 63.5463
    "Zn" : 65.4094,     # 30 Zn zinc 65.4094
    "Ga" : 69.7231,     # 31 Ga gallium 69.7231
    "Ge" : 72.641,      # 32 Ge germanium 72.641
    "As" : 74.921602,   # 33 As arsenic 74.921602
    "Se" : 78.963,      # 34 Se selenium 78.963
    "Br" : 79.9041,     # 35 Br bromine 79.9041
    "Kr" : 83.7982,     # 36 Kr krypton 83.7982
    "Rb" : 85.46783,    # 37 Rb rubidium 85.46783
    "Sr" : 87.621,      # 38 Sr strontium 87.621
    "Y"  : 88.905852,   # 39 Y yttrium 88.905852
    "Zr" : 91.2242,     # 40 Zr zirconium 91.2242
    "Nb" : 92.906382,   # 41 Nb niobium 92.906382
    "Mo" : 95.942,      # 42 Mo molybdenum 95.942
    "Tc" : 98.0,        # 43 Tc technetium 98
    "Ru" : 101.072,     # 44 Ru ruthenium 101.072
    "Rh" : 102.905502,  # 45 Rh rhodium 102.905502
    "Pd" : 106.421,     # 46 Pd palladium 106.421
    "Ag" : 107.86822,   # 47 Ag silver 107.86822
    "Cd" : 112.4118,    # 48 Cd cadmium 112.4118
    "In" : 114.8183,    # 49 In indium 114.8183
    "Sn" : 118.7107,    # 50 Sn tin 118.7107
    "Sb" : 121.7601,    # 51 Sb antimony 121.7601
    "Te" : 127.603,     # 52 Te tellurium 127.603
    "I"  : 126.904473,  # 53 I iodine 126.904473
    "Xe" : 131.2936,    # 54 Xe xenon 131.2936
    "Cs" : 132.905452,  # 55 Cs caesium 132.905452
    "Ba" : 137.3277,    # 56 Ba barium 137.3277
    "La" : 138.90552,   # 57 La lanthanum 138.90552
    "Ce" : 140.1161,    # 58 Ce cerium 140.1161
    "Pr" : 140.907652,  # 59 Pr praseodymium 140.907652
    "Nd" : 144.243,     # 60 Nd neodymium 144.243
    "Pm" : 145.0,       # 61 Pm promethium 145
    "Sm" : 150.363,     # 62 Sm samarium 150.363
    "Eu" : 151.9641,    # 63 Eu europium 151.9641
    "Gd" : 157.253,     # 64 Gd gadolinium 157.253
    "Tb" : 158.925342,  # 65 Tb terbium 158.925342
    "Dy" : 162.5001,    # 66 Dy dysprosium 162.5001
    "Ho" : 164.930322,  # 67 Ho holmium 164.930322
    "Er" : 167.2593,    # 68 Er erbium 167.2593
    "Tm" : 168.934212,  # 69 Tm thulium 168.934212
    "Yb" : 173.043,     # 70 Yb ytterbium 173.043
    "Lu" : 174.9671,    # 71 Lu lutetium 174.9671
    "Hf" : 178.492,     # 72 Hf hafnium 178.492
    "Ta" : 180.94791,   # 73 Ta tantalum 180.94791
    "W"  : 183.841,     # 74 W tungsten 183.841
    "Re" : 186.2071,    # 75 Re rhenium 186.2071
    "Os" : 190.233,     # 76 Os osmium 190.233
    "Ir" : 192.2173,    # 77 Ir iridium 192.2173
    "Pt" : 195.0782,    # 78 Pt platinum 195.0782
    "Au" : 196.966552,  # 79 Au gold 196.966552
    "Hg" : 200.592,     # 80 Hg mercury 200.592
    "Tl" : 204.38332,   # 81 Tl thallium 204.38332
    "Pb" : 207.21,      # 82 Pb lead 207.21
    "Bi" : 208.980382,  # 83 Bi bismuth 208.980382
    "Po" : 209.0,       # 84 Po polonium 209
    "At" : 210.0,       # 85 At astatine 210
    "Rn" : 222.0,       # 86 Rn radon 222
    "Fr" : 223.0,       # 87 Fr francium 223
    "Ra" : 226.0,       # 88 Ra radium 226
    "Ac" : 227.0,       # 89 Ac actinium 227
    "Th" : 232.03811,   # 90 Th thorium 232.03811
    "Pa" : 231.035882,  # 91 Pa protactinium 231.035882
    "U"  : 238.028913,  # 92 U uranium 238.028913
    "Np" : 237.0,       # 93 Np neptunium 237
    "Pu" : 244.0,       # 94 Pu plutonium 244
    "Am" : 243.0,       # 95 Am americium 243
    "Cm" : 247.0,       # 96 Cm curium 247
    "Bk" : 247.0,       # 97 Bk berkelium 247
    "Cf" : 251.0,       # 98 Cf californium 251
    "Es" : 252.0,       # 99 Es einsteinium 252
    "Fm" : 257.0,       # 100 Fm fermium 257
    "Md" : 258.0,       # 101 Md mendelevium 258
    "No" : 259.0,       # 102 No nobelium 259
    "Lr" : 262.0,       # 103 Lr lawrencium 262
    "Rf" : 261.0,       # 104 Rf rutherfordium 261
    "Db" : 262.0,       # 105 Db dubnium 262
    "Sg" : 266.0,       # 106 Sg seaborgium 266
    "Bh" : 264.0,       # 107 Bh bohrium 264
    "Hs" : 277.0,       # 108 Hs hassium 277
    "Mt" : 268.0,       # 109 Mt meitnerium 268
    "Ds" : 281.0,       # 110 Ds darmstadtium 281
    "Rg" : 272.0,       # 111 Rg roentgenium 272
}

# End of Constants


##############################################################################
class P_xcfg(StructureParser):
    """Parser for AtomEye extended CFG format.

    cluster_boundary -- width of boundary around corners of non-periodic
                        cluster to avoid PBC effects in atomeye
    """

    cluster_boundary = 2

    def __init__(self):
        StructureParser.__init__(self)
        self.format = "xcfg"
        return

    def parseLines(self, lines):
        """Parse list of lines in PDB format.

        Return Structure object or raise StructureFormatError.
        """
        xcfg_Number_of_particles = None
        xcfg_A = None
        xcfg_H0 = numpy.zeros((3,3), dtype=float)
        xcfg_H0_set = numpy.zeros((3,3), dtype=bool)
        xcfg_NO_VELOCITY = False
        xcfg_entry_count = None
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
                        emsg = ("%d: first line must " +
                                "contain 'Number of particles ='") % p_nl
                        raise StructureFormatError(emsg)
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
            if numpy.any(xcfg_H0_set == False):
                emsg = "H0 tensor is not properly defined"
                raise StructureFormatError(emsg)
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
                emsg = ("%d: auxiliary fields " +
                        "not consistent with entry_count") % p_nl
                raise StructureFormatError(emsg)
            # define proper lattice
            stru.lattice.setLatBase(xcfg_H0)
            # build p_assign_atom function to assign entries to proper fields
            p_exprs = [ "a.xyz[0]=fields[0]",
                        "a.xyz[1]=fields[1]",
                        "a.xyz[2]=fields[2]" ]
            if not xcfg_NO_VELOCITY:
                p_exprs += [  "a.v=numpy.zeros(3, dtype=float)",
                              "a.v[0]=fields[3]",
                              "a.v[1]=fields[4]",
                              "a.v[2]=fields[5]" ]
            for idx in sorted_aux_keys:
                prop = p_auxiliary[idx]
                col = idx + 6 - 3*xcfg_NO_VELOCITY
                if prop == "Uiso":
                    p_exprs.append("a.Uisoequiv=fields[%d]" % col)
                elif re.match(r"^U\d\d$", prop) \
                and 1<=int(prop[1])<=3 and 1<=int(prop[2])<=3 :
                    p_exprs.append("a.anisotropy=True")
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
                # ignore atom mass
                if len(words) == 1 and isfloat(words[0]):
                    continue
                # parse element allowing empty symbol
                elif len(words) <= 1:
                    w = line.strip()
                    p_element = w[:1].upper() + w[1:].lower()
                elif len(words) == xcfg_entry_count and p_element is not None:
                    fields = [ float(w) for w in words ]
                    stru.addNewAtom(p_element, fields[:3])
                    a = stru.getLastAtom()
                    a.xyz *= xcfg_A
                    p_assign_atom(a, fields)
                else:
                    emsg = "%d: invalid record" % p_nl
                    raise StructureFormatError(emsg)
            if len(stru) != p_natoms:
                emsg = "expected %d atoms, read %d" % (p_natoms, len(stru))
                raise StructureFormatError(emsg)
        except (ValueError, IndexError):
            emsg = "%d: file is not in XCFG format" % p_nl
            exc_type, exc_value, exc_traceback = sys.exc_info()
            raise StructureFormatError, emsg, exc_traceback
        return stru
    # End of parseLines

    def toLines(self, stru):
        """Convert Structure stru to a list of lines in XCFG atomeye format.

        Return list of strings.
        """
        if len(stru) == 0:
            emsg = "cannot convert empty structure to XCFG format"
            raise StructureFormatError(emsg)
        lines = []
        lines.append( "Number of particles = %i" % len(stru) )
        # figure out length unit A
        allxyz = numpy.array([a.xyz for a in stru])
        lo_xyz = allxyz.min(axis=0)
        hi_xyz = allxyz.max(axis=0)
        max_range_xyz = (hi_xyz-lo_xyz).max()
        if numpy.allclose(stru.lattice.abcABG(), (1, 1, 1, 90, 90, 90)):
            max_range_xyz += self.cluster_boundary
        # range of CFG coordinates must be less than 1
        p_A = numpy.ceil(max_range_xyz + 1.0e-13)
        # atomeye draws rubbish when boxsize is less than 3.5
        hi_ucvect = max([numpy.sqrt(numpy.dot(v,v)) for v in stru.lattice.base])
        if hi_ucvect*p_A < 3.5:
            p_A = numpy.ceil(3.5 / hi_ucvect)
        lines.append( "A = %.8g Angstrom" % p_A )
        # how much do we need to shift the coordinates?
        p_dxyz = numpy.zeros(3, dtype=float)
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
            if p_allUzero and numpy.any(a.U != 0.0):
                p_allUzero = False
            if not numpy.all(a.U == a.U[0,0]*numpy.identity(3)):
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
            allU = numpy.array([a.U for a in stru])
            if numpy.any(allU[:,0,1] != 0.0):
                p_auxiliaries.append( ('U12', 'a.U[0,1]') )
            if numpy.any(allU[:,0,2] != 0.0):
                p_auxiliaries.append( ('U13', 'a.U[0,2]') )
            if numpy.any(allU[:,1,2] != 0.0):
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

# End of class P_xcfg

# Routines

def getParser():
    return P_xcfg()

# End of file

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

"""Parser for PDFFit file format"""

__id__ = "$Id$"

import sys
import numpy as num

from Structure.PDFFitStructure import PDFFitStructure
from Structure.Lattice import Lattice
from Structure.Atom import Atom
from StructureParser import StructureParser
from Structure.exceptions import InvalidStructureFormat

class Parser(StructureParser):
    """Parser --> StructureParser subclass for PDFFit format"""

    def __init__(self):
        self.format = "pdffit"
        return

    def parseLines(self, lines):
        """Parse list of lines in PDFFit format.

        Return Structure object or raise InvalidStructureFormat.
        """
        p_nl = 0
        rlist = []
        try:
            stru = PDFFitStructure()
            cell_line_read = False
            stop = len(lines)
            while stop>0 and lines[stop-1].strip() == "":
                stop -= 1
            ilines = iter(lines[:stop])
            # read header of PDFFit file
            for l in ilines:
                p_nl += 1
                words = l.split()
                if len(words) == 0 or words[0][0] == '#':
                    continue
                elif words[0] == 'title':
                    stru.title = l.lstrip()[5:].strip()
                elif words[0] == 'scale':
                    stru.pdffit['scale'] = float(words[1])
                elif words[0] == 'sharp':
                    l1 = l.replace(',', ' ')
                    sharp_pars = [ float(w) for w in l1.split()[1:] ]
                    if len(sharp_pars) < 4:
                        stru.pdffit['delta'] = sharp_pars[0]
                        stru.pdffit['srat'] = sharp_pars[1]
                        stru.pdffit['rcut'] = sharp_pars[2]
                    else:
                        stru.pdffit['delta'] = sharp_pars[0]
                        stru.pdffit['gamma'] = sharp_pars[1]
                        stru.pdffit['srat'] = sharp_pars[2]
                        stru.pdffit['rcut'] = sharp_pars[3]
                elif words[0] == 'spcgr':
                    stru.pdffit['spcgr'] = ''.join(words[1:])
                elif words[0] == 'cell':
                    cell_line_read = True
                    l1 = l.replace(',', ' ')
                    latpars = [ float(w) for w in l1.split()[1:7] ]
                    stru.lattice = Lattice(*latpars)
                elif words[0] == 'dcell':
                    l1 = l.replace(',', ' ')
                    stru.pdffit['dcell'] = [ float(w) for w in l1.split()[1:7] ]
                elif words[0] == 'ncell':
                    l1 = l.replace(',', ' ')
                    stru.pdffit['ncell'] = [ int(w) for w in l1.split()[1:5] ]
                elif words[0] == 'format':
                    if words[1] != 'pdffit':
                        raise InvalidStructureFormat, \
                                "%d: file is not in PDFFit format" % p_nl
                elif words[0] == 'atoms' and cell_line_read:
                    break
                else:
                    raise InvalidStructureFormat, \
                            "%d: file is not in PDFFit format" % p_nl
            # header was successfully read, we can create Structure instance
            p_natoms = reduce(lambda x,y : x*y, stru.pdffit['ncell'])
            # we are now inside data block
            for l in ilines:
                p_nl += 1
                wl1 = l.split()
                element = wl1[0][0].upper() + wl1[0][1:].lower()
                xyz = [ float(w) for w in wl1[1:4] ]
                a = Atom(element, xyz=xyz, occupancy=float(wl1[4]))
                p_nl += 1
                wl2 = ilines.next().split()
                a.sigxyz = [ float(w) for w in wl2[0:3] ]
                a.sigo = float(wl2[3])
                p_nl += 1
                wl3 = ilines.next().split()
                p_nl += 1
                wl4 = ilines.next().split()
                p_nl += 1
                wl5 = ilines.next().split()
                p_nl += 1
                wl6 = ilines.next().split()
                a.sigU = num.zeros((3,3), dtype=float)
                for i in range(3):
                    a.U[i][i] = float(wl3[i])
                    a.sigU[i][i] = float(wl4[i])
                a.U[0][1] = a.U[1][0] = float(wl5[0])
                a.U[0][2] = a.U[2][0] = float(wl5[1])
                a.U[1][2] = a.U[2][1] = float(wl5[2])
                a.sigU[0][1] = a.sigU[1][0] = float(wl6[0])
                a.sigU[0][2] = a.sigU[2][0] = float(wl6[1])
                a.sigU[1][2] = a.sigU[2][1] = float(wl6[2])
                stru.append(a)
            if len(stru) != p_natoms:
                raise InvalidStructureFormat, \
                        "expected %d atoms, read %d" % (p_natoms, len(stru))
            if stru.pdffit['ncell'][:3] != [1,1,1]:
                superlatpars = [ latpars[i]*stru.pdffit['ncell'][i]
                                 for i in range(3) ] + latpars[3:]
                superlattice = Lattice(*superlatpars)
                stru.placeInLattice(superlattice)
                stru.pdffit['ncell'] = [1, 1, 1, p_natoms]
        except (ValueError, IndexError):
            exc_type, exc_value, exc_traceback = sys.exc_info()
            raise InvalidStructureFormat, \
                    "%d: file is not in PDFFit format" % p_nl, exc_traceback
        return stru
    # End of parseLines

    def toLines(self, stru):
        """Convert Structure stru to a list of lines in PDFFit format.

        Return list of strings.
        """
        # first, convert stru to PDFFitStructure
        if not isinstance(stru, PDFFitStructure):
            pfstru = PDFFitStructure()
            pfstru.__dict__.update(stru.__dict__)
            pfstru[:] = stru[:]
            stru = pfstru
        lines = []
        # default values of standard deviations
        d_sigxyz = num.zeros(3, dtype=float)
        d_sigo = 0.0
        d_sigU = num.zeros((3,3), dtype=float)
        # here we can start
        l = "title  " + stru.title
        lines.append( l.strip() )
        lines.append( "format pdffit" )
        lines.append( "scale  %9.6f" % stru.pdffit["scale"] )
        lines.append( "sharp  %9.6f, %9.6f, %9.6f, %9.6f" % (
            stru.pdffit["delta"],
            stru.pdffit["gamma"],
            stru.pdffit["srat"],
            stru.pdffit["rcut"]) )
        lines.append( "spcgr   " + stru.pdffit["spcgr"] )
        lat = stru.lattice
        lines.append( "cell   %9.6f, %9.6f, %9.6f, %9.6f, %9.6f, %9.6f" % (
            lat.a, lat.b, lat.c, lat.alpha, lat.beta, lat.gamma) )
        lines.append( "dcell  %9.6f, %9.6f, %9.6f, %9.6f, %9.6f, %9.6f" %
            tuple(stru.pdffit["dcell"]) )
        lines.append( "ncell  %9i, %9i, %9i, %9i" % (1, 1, 1, len(stru)) )
        lines.append( "atoms" )
        for a in stru:
            ad = a.__dict__
            lines.append( "%-4s%18.8f%18.8f%18.8f%13.4f" % (
                a.element.upper(), a.xyz[0], a.xyz[1], a.xyz[2], a.occupancy) )
            sigmas = num.concatenate(
                ( ad.get("sigxyz", d_sigxyz),  [ad.get("sigo", d_sigo)] )  )
            lines.append( "    %18.8f%18.8f%18.8f%13.4f" % tuple(sigmas) )
            sigU = ad.get("sigU", d_sigU)
            Uii = ( a.U[0][0], a.U[1][1], a.U[2][2] )
            Uij = ( a.U[0][1], a.U[0][2], a.U[1][2] )
            sigUii = ( sigU[0][0], sigU[1][1], sigU[2][2] )
            sigUij = ( sigU[0][1], sigU[0][2], sigU[1][2] )
            lines.append( "    %18.8f%18.8f%18.8f" % Uii )
            lines.append( "    %18.8f%18.8f%18.8f" % sigUii )
            lines.append( "    %18.8f%18.8f%18.8f" % Uij )
            lines.append( "    %18.8f%18.8f%18.8f" % sigUij )
        return lines
    # End of toLines

# End of Parser

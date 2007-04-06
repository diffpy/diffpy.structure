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

"""This module defines class Structure.
"""

__id__ = "$Id$"

import copy
import math
import numpy
import numpy.linalg as numalg
from Lattice import Lattice
from Atom import Atom

##############################################################################
class Structure(list):
    """Structure --> group of atoms

    Structure class is inherited from Python list.  It contains
    a list of Atom instances.

    Data members:
        title   -- structure description
        lattice -- coordinate system (instance of Lattice)
        fileformat -- last format used for reading or writing of file
    """

    def __init__(self, atoms=[], lattice=None, title=""):
        """define group of atoms in a specified lattice.

        atoms   -- list of Atom instances to be included in this Structure.
                   When atoms argument is an existing Structure instance,
                   the new Structure is its deep copy.
        lattice -- instance of Lattice that defines coordinate system
        title   -- string description of the structure

        Structure(stru)     create a copy of Structure instance stru.

        Because Structure is inherited from a list it can use list expansions,
        for example:
            oxygen_atoms = [ for a in stru if a.element == "O" ]
            oxygen_stru = Structure(oxygen_atoms, lattice=stru.lattice)
        """
        self.extend(atoms)
        self.title = ""
        self.lattice = None
        self.fileformat = None
        if isinstance(atoms, Structure):
            stru = atoms
            self.extend([Atom(a) for a in stru])
            self.lattice = Lattice(stru.lattice)
            self.title = stru.title
        if lattice is None:
            if not self.lattice:    self.lattice = Lattice()
        elif not isinstance(lattice, Lattice):
            raise TypeError, "expected instance of Lattice"
        else:
            self.lattice = lattice
        return

    def __str__(self):
        """simple string representation"""
        s_lattice = "lattice=%s" % self.lattice
        s_atoms = '\n'.join([str(a) for a in self])
        return s_lattice + '\n' + s_atoms

    def dist(self, a0, a1):
        """distance of 2 atoms"""
        return self.lattice.dist(a0.xyz, a1.xyz)

    def angle(self, a0, a1, a2):
        """angle at atom a1 in degrees"""
        u10 = [ (d[1]-d[0]) for d in zip(a0.xyz, a1.xyz) ]
        u12 = [ (d[1]-d[0]) for d in zip(a2.xyz, a1.xyz) ]
        return self.lattice.angle(u10, u12)

    def cartesian(self, a):
        """return cartesian coordinates of atom a"""
        vc = self.lattice.cartesian(a.xyz)
        return vc

    def msdLat(self, a, vl):
        """mean square displacement of an atom along lattice vector

        a  -- atom instance
        vl -- vector in lattice coordinates

        return mean square displacement
        """
        vln = numpy.array(vl, dtype=float)/self.lattice.norm(vl)
        G = self.lattice.metrics
        rhs = numpy.array([ G[0]*self.lattice.ar,
                          G[1]*self.lattice.br,
                          G[2]*self.lattice.cr ], dtype=float)
        rhs = numpy.dot(rhs, vln)
        msd = numpy.dot(rhs, numpy.dot(a.U, rhs))
        return msd

    def msdCart(self, a, vc):
        """mean square displacement of an atom along cartesian vector

        a  -- atom instance
        vc -- vector in absolute cartesian coordinates

        return mean square displacement
        """
        vcn = numpy.array(vc, dtype=float)
        vcn /= numpy.sqrt(numpy.sum(vcn**2))
        F1 = self.lattice.normbase
        Uc = numpy.dot(numpy.transpose(F1), numpy.dot(a.U, F1))
        msd = numpy.dot(vcn, numpy.dot(Uc, vcn))
        return msd

    def placeInLattice(self, new_lattice):
        """place structure into new_lattice coordinate system

        sets lattice to new_lattice and recalculate fractional coordinates
        of all atoms so their absolute positions remain the same

        return self
        """
        Tx = numpy.dot(self.lattice.base, new_lattice.recbase)
        Tu = numpy.dot(self.lattice.normbase, new_lattice.recnormbase)
        for a in self:
            a.xyz = numpy.dot(a.xyz, Tx)
            a.U = numpy.dot(numpy.transpose(Tu), numpy.dot(a.U, Tu))
        self.lattice = new_lattice
        return self

    def read(self, filename, format='auto'):
        """load structure from file, original data get lost

        filename -- file to be loaded
        format   -- all structure formats are defined in Parsers submodule,
                    when format == 'auto' all Parsers are tried one by one

        return self
        """
        self.readStr(open(filename,'r').read(), format)
        if not self.title:
            import os.path
            self.title = os.path.basename(filename)
            self.title = os.path.splitext(self.title)[0]
        return self

    def readStr(self, s, format='auto'):
        """load structure from a string, original data get lost

        s        -- string with structure definition
        format   -- all structure formats are defined in Parsers submodule,
                    when format == 'auto' all Parsers are tried one by one

        return self
        """
        from Parsers import parse
        new_structure = parse(s, format)
        self.__dict__.update(new_structure.__dict__)
        self[:] = new_structure[:]
        return self

    def write(self, filename, format):
        """save structure to file in the specified format

        Note: available structure formats can be obtained by:
            from Parsers import formats
        """
        s = self.writeStr(format)
        file = open(filename,'w')
        file.write(s)
        file.close()
        return

    def writeStr(self, format):
        """return string representation of the structure in specified format

        Note: available structure formats can be obtained by:
            from Parsers import formats
        """
        from Parsers import tostring
        return tostring(self, format)

# End of class Structure

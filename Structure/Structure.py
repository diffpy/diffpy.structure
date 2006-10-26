"""This module defines class Structure.
"""

__id__ = "$Id$"

import copy
import math
import numpy as num
import numpy.linalg as numalg
from Lattice import Lattice
from Atom import Atom

##############################################################################
class Structure(list):
    """Structure --> group of atoms

    Structure class is inherited from Python list

    Data members:
        title   -- structure description
        lattice -- coordinate system (instance of Lattice)
    """

    def __init__(self, atoms=[], lattice=None, title=""):
        """define group of atoms in a specified lattice

        atoms   -- list of Atom instances, when atoms is Structure,
                   create a deep copy of atoms
        lattice -- instance of Lattice that defines coordinate system
        title   -- string title

        because Structure is inherited from a list it can use list expansions,
        for example:
            oxygen_atoms = [ for a in stru if a.element == "O" ]
            oxygen_stru = Structure(oxygen_atoms, lattice=stru.lattice)
        """
        self.extend(atoms)
        self.title = title
        if isinstance(atoms, Structure):
            for attribute, value in atoms.__dict__.iteritems():
                setattr(self, attribute, copy.deepcopy(value))
        elif lattice is None:
            self.lattice = Lattice()
        elif not isinstance(lattice, Lattice):
            raise RuntimeError, "expected instance of Lattice"
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
        vln = num.array(vl, dtype=float)/self.lattice.norm(vl)
        G = self.lattice.metrics
        rhs = num.array([ G[0]*self.lattice.ar,
                          G[1]*self.lattice.br,
                          G[2]*self.lattice.cr ], dtype=float)
        rhs = num.dot(rhs, vln)
        msd = num.dot(rhs, num.dot(a.U, rhs))
        return msd

    def msdCart(self, a, vc):
        """mean square displacement of an atom along cartesian vector

        a  -- atom instance
        vc -- vector in absolute cartesian coordinates

        return mean square displacement
        """
        vcn = num.array(vc, dtype=float)
        vcn /= num.sqrt(num.sum(vcn**2))
        F1 = self.lattice.normbase
        Uc = num.dot(num.transpose(F1), num.dot(a.U, F1))
        msd = num.dot(vcn, num.dot(Uc, vcn))
        return msd

    def placeInLattice(self, new_lattice):
        """place structure into new_lattice coordinate system

        sets lattice to new_lattice and recalculate fractional coordinates
        of all atoms so their absolute positions remain the same

        return self
        """
        Tx = num.dot(self.lattice.base, new_lattice.recbase)
        Tu = num.dot(self.lattice.normbase, new_lattice.recnormbase)
        for a in self:
            a.xyz = num.dot(a.xyz, Tx)
            a.U = num.dot(num.transpose(Tu), num.dot(a.U, Tu))
        self.lattice = new_lattice
        return self

    def read(self, filename, format='auto'):
        """load structure from file, original data get lost

        filename -- self explanatory
        format   -- all structure formats are defined in Parsers submodulede,
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
        format   -- all structure formats are defined in Parsers submodulede,
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

# End of Structure

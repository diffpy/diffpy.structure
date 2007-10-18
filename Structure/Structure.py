########################################################################
#
# Structure         by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2007 trustees of the Michigan State University.
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
    a list of Atom instances.  Structure overloads setitem and setslice
    methods so that the lattice attribute of atoms get set to lattice.

    Data members:
        title   -- structure description
        lattice -- coordinate system (instance of Lattice)
    """

    def __init__(self, atoms=[], lattice=None, title=""):
        """define group of atoms in a specified lattice.

        atoms   -- list of Atom instances to be included in this Structure.
                   When atoms argument is an existing Structure instance,
                   the new Structure is its deep copy.
        lattice -- instance of Lattice defining coordinate systems, property.
        title   -- string description of the structure

        Structure(stru)     create a copy of Structure instance stru.

        Because Structure is inherited from a list it can use list expansions,
        for example:
            oxygen_atoms = [ for a in stru if a.element == "O" ]
            oxygen_stru = Structure(oxygen_atoms, lattice=stru.lattice)
        """
        self.title = ""
        self._lattice = None
        if isinstance(atoms, Structure):
            # copy lattice and title
            stru = atoms
            self.lattice = Lattice(stru.lattice)
            self.title = stru.title
        # override from lattice argument
        if lattice is None:
            if not self.lattice:    self.lattice = Lattice()
        elif not isinstance(lattice, Lattice):
            raise TypeError, "expected instance of Lattice"
        else:
            self.lattice = lattice
        # override from title argument
        if title:
            self.title = title
        # finally assign list of atoms to self
        self[:] = atoms
        return

    def __str__(self):
        """simple string representation"""
        s_lattice = "lattice=%s" % self.lattice
        s_atoms = '\n'.join([str(a) for a in self])
        return s_lattice + '\n' + s_atoms

    def addNewAtom(self, *args, **kwargs):
        """Add new Atom instance to the end of this Structure.

        All arguments are forwarded to Atom constructor.

        No return value.
        """
        kwargs['lattice'] = self.lattice
        a = Atom(*args, **kwargs)
        list.append(self, a)
        return

    def getLastAtom(self):
        """Return Reference to the last Atom in this structure.
        """
        last_atom = self[-1]
        return last_atom

    def dist(self, a0, a1):
        """distance of 2 atoms"""
        return self.lattice.dist(a0.xyz, a1.xyz)

    def angle(self, a0, a1, a2):
        """angle at atom a1 in degrees"""
        u10 = [ (d[1]-d[0]) for d in zip(a0.xyz, a1.xyz) ]
        u12 = [ (d[1]-d[0]) for d in zip(a2.xyz, a1.xyz) ]
        return self.lattice.angle(u10, u12)

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
            if a.anisotropy:
                a.U = numpy.dot(numpy.transpose(Tu), numpy.dot(a.U, Tu))
        self.lattice = new_lattice
        return self

    def read(self, filename, format='auto'):
        """Load structure from a file, any original data become lost.

        filename -- file to be loaded
        format   -- all structure formats are defined in Parsers submodule,
                    when format == 'auto' all Parsers are tried one by one

        Return instance of data Parser used to process file.  This
        can be inspected for information related to particular format.
        """
        from Parsers import getParser
        p = getParser(format)
        new_structure = p.parseFile(filename)
        # reinitialize data after successful parsing
        # avoid calling __init__ from a derived class
        Structure.__init__(self)
        if new_structure is not None:
            self.__dict__.update(new_structure.__dict__)
            self[:] = new_structure
        if not self.title:
            import os.path
            tailname = os.path.basename(filename)
            tailbase = os.path.splitext(tailname)[0]
            self.title = tailbase
        return p

    def readStr(self, s, format='auto'):
        """Load structure from a string, any original data become lost.

        s        -- string with structure definition
        format   -- all structure formats are defined in Parsers submodule,
                    when format == 'auto' all Parsers are tried one by one

        Return instance of data Parser used to process input string.  This
        can be inspected for information related to particular format.
        """
        from Parsers import getParser
        p = getParser(format)
        new_structure = p.parse(s)
        # reinitialize data after successful parsing
        # avoid calling __init__ from a derived class
        Structure.__init__(self)
        if new_structure is not None:
            self.__dict__.update(new_structure.__dict__)
            self[:] = new_structure
        return p

    def write(self, filename, format):
        """save structure to file in the specified format

        Note: available structure formats can be obtained by:
            from Parsers import formats
        """
        from Parsers import getParser
        p = getParser(format)
        p.filename = filename
        s = p.tostring(self)
        f = open(filename, 'wb')
        f.write(s)
        f.close()
        return

    def writeStr(self, format):
        """return string representation of the structure in specified format

        Note: available structure formats can be obtained by:
            from Parsers import formats
        """
        from Parsers import getParser
        p = getParser(format)
        s = p.tostring(self)
        return s

    ##############################################################################
    # overloaded list methods
    ##############################################################################

    def append(self, a, copy=True):
        """Append atom to a structure and update its lattice attribute.

        a    -- instance of Atom
        copy -- flag for appending a copy of a.  
                When False, append a and update a.owner.

        No return value.
        """
        adup = copy and Atom(a) or a
        adup.lattice = self.lattice
        list.append(self, adup)
        return

    def insert(self, idx, a, copy=True):
        """Insert atom a before position idx in this Structure.

        idx  -- position in atom list
        a    -- instance of Atom
        copy -- flag for inserting a copy of a.  
                When False, append a and update a.lattice.

        No return value.
        """
        adup = copy and Atom(a) or a
        adup.lattice = self.lattice
        list.insert(idx, adup)
        return

    def extend(self, atoms, copy=True):
        """Extend Structure by appending copies from a list of atoms.

        atoms -- list of Atom instances
        copy  -- flag for extending with copies of Atom instances.
                 When False extend with atoms and update their lattice
                 attributes.

        No return value.
        """
        if copy:    adups = [Atom(a) for a in atoms]
        else:       adups = atoms
        for a in adups: a.lattice = self.lattice
        list.extend(self, adups)
        return

    def __setitem__(self, idx, a, copy=True):
        """Set idx-th atom to a.

        idx  -- index of atom in this Structure
        a    -- instance of Atom
        copy -- flag for setting to a copy of a.  
                When False, set to a and update a.lattice.

        No return value.
        """
        adup = copy and Atom(a) or a
        adup.lattice = self.lattice
        list.__setitem__(self, idx, adup)
        return

    def __setslice__(self, lo, hi, atoms):
        """Set Structure slice from lo to hi-1 to the sequence of atoms.

        lo    -- low index for the slice
        hi    -- high index of the slice
        atoms -- sequence of Atom instances
        copy  -- flag for using copies of Atom instances.  When False, set
                 to existing instances and update their lattice attributes.

        No return value.
        """
        if copy:    adups = [Atom(a) for a in atoms]
        else:       adups = atoms
        for a in adups: a.lattice = self.lattice
        list.__setslice__(self, lo, hi, adups)
        return


    ####################################################################
    # property handlers
    ####################################################################

    # lattice

    def _get_lattice(self):
        return self._lattice

    def _set_lattice(self, value):
        for a in self:  a.lattice = value
        self._lattice = value
        return

    lattice = property(_get_lattice, _set_lattice, doc =
        "Coordinate system for this Structure.")


# End of class Structure

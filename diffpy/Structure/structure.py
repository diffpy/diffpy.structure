##############################################################################
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
##############################################################################

"""This module defines class Structure.
"""

__id__ = "$Id$"

import numpy
import copy
from diffpy.Structure import Lattice, Atom
from diffpy.Structure.utils import _linkAtomAttribute

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

    def __init__(self, atoms=[], lattice=None, title="", filename=None,
            format=None):
        """define group of atoms in a specified lattice.

        atoms    -- list of Atom instances to be included in this Structure.
                    When atoms argument is an existing Structure instance,
                    the new Structure is its copy.
        lattice  -- instance of Lattice defining coordinate systems, property.
        title    -- string description of the structure
        filename -- optional, name of a file to load the structure from.
                    Overrides atoms argument when specified.
        format   -- optional structure format of the loaded filename.  By default
                    all structure formats are tried one by one.  Ignored when
                    filename has not been specified.

        Structure(stru)     create a copy of Structure instance stru.

        Because Structure is inherited from a list it can use list expansions,
        for example:
            oxygen_atoms = [ for a in stru if a.element == "O" ]
            oxygen_stru = Structure(oxygen_atoms, lattice=stru.lattice)
        """
        self.title = ""
        self._lattice = None
        self._labels = {}
        self._labels_cached = False
        if isinstance(atoms, Structure):
            atoms.__copy__(target=self)
        # override from lattice argument
        if lattice is None:
            if not self.lattice:    self.lattice = Lattice()
        elif not isinstance(lattice, Lattice):
            emsg = "expected instance of Lattice"
            raise TypeError(emsg)
        else:
            self.lattice = lattice
        # override from title argument
        if title:
            self.title = title
        # finally check if data should be loaded from file
        if filename is not None:
            readkwargs = {}
            if format is not None:  readkwargs['format'] = format
            self.read(filename, **readkwargs)
        # otherwise assign list of atoms to self unless already done by copy
        elif len(self) != len(atoms):
            self[:] = atoms
        return


    def __copy__(self, target=None, shallow=False):
        '''Create a deep copy of this instance.

        target   -- optional target instance for copying, useful for
                    copying a derived class.  Defaults to new instance
                    of the same type as self.
        shallow  -- flag for creating a shallow copy

        Return a duplicate instance of this object.
        '''
        if target is None:
            target = Structure()
        elif target is self:
            return target
        # create a shallow copy of all source attributes
        target.__dict__.update(self.__dict__)
        if not shallow:
            # make a deep copy of the lattice and pdffit dictionary
            target.lattice = Lattice(self.lattice)
            if hasattr(self, 'pdffit'):
                target.pdffit = copy.deepcopy(self.pdffit)
        # copy all atoms to the target
        target.__setslice__(0, len(target), self, copy=not shallow)
        return target


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
        self._uncache('labels')
        return


    def getLastAtom(self):
        """Return Reference to the last Atom in this structure.
        """
        last_atom = self[-1]
        return last_atom


    def getAtom(self, id):
        """Reference to internal Atom specified by the identifier.

        id  -- zero based index or a string label formatted as
               "%(element)s%(order)i", for example "Na1", "Cl1"

        Return Atom instance.
        Raise ValueError for invalid id.

        See also getLabels().
        """
        try:
            if type(id) is int:
                rv = self[id]
            else:
                if not self._labels_cached or id not in self._labels:
                    self._update_labels()
                rv = self._labels[id]
        except (IndexError, KeyError):
            emsg = "Invalid atom identifier %r." % id
            raise ValueError(emsg)
        return rv


    def getLabels(self):
        """List of unique string labels for all atoms in this structure.

        Return a list.
        """
        elnum = {}
        labels = []
        for a in self:
            elnum[a.element] = elnum.get(a.element, 0) + 1
            alabel = a.element + str(elnum[a.element])
            labels.append(alabel)
        return labels


    def distance(self, id0, id1):
        """Distance between 2 atoms, no periodic boundary conditions.

        id0  -- zero based index of the first atom or a string label
                such as "Na1"
        id1  -- zero based index or string label of the second atom.

        Return float.
        Raise ValueError for invalid arguments.
        """
        a0 = self.getAtom(id0)
        a1 = self.getAtom(id1)
        return self.lattice.dist(a0.xyz, a1.xyz)


    def angle(self, id0, id1, id2):
        """The bond angle at the second of three atoms in degrees.

        id0  -- zero based index of the first atom or a string label
                such as "Na1"
        id1  -- index or string label for the second atom,
                where the angle is formed
        id2  -- index or string label for the third atom

        Return float.
        Raise ValueError for invalid arguments.
        """
        a0 = self.getAtom(id0)
        a1 = self.getAtom(id1)
        a2 = self.getAtom(id2)
        u10 = a0.xyz - a1.xyz
        u12 = a2.xyz - a1.xyz
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
        import diffpy.Structure
        import diffpy.Structure.Parsers
        getParser = diffpy.Structure.Parsers.getParser
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
        from diffpy.Structure.Parsers import getParser
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
        """Save structure to file in the specified format

        No return value.

        Note: available structure formats can be obtained by:
            from Parsers import formats
        """
        from diffpy.Structure.Parsers import getParser
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
        from diffpy.Structure.Parsers import getParser
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
                When False, append a and update a.lattice.

        No return value.
        """
        self._uncache('labels')
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
        self._uncache('labels')
        adup = copy and Atom(a) or a
        adup.lattice = self.lattice
        list.insert(self, idx, adup)
        return


    def extend(self, atoms, copy=True):
        """Extend Structure by appending copies from a list of atoms.

        atoms -- list of Atom instances
        copy  -- flag for extending with copies of Atom instances.
                 When False extend with atoms and update their lattice
                 attributes.

        No return value.
        """
        self._uncache('labels')
        if copy:    adups = [Atom(a) for a in atoms]
        else:       adups = atoms
        for a in adups: a.lattice = self.lattice
        list.extend(self, adups)
        return


    def __getitem__(self, idx):
        try:
            return list.__getitem__(self, idx)
        except TypeError:
            pass
        indices = numpy.arange(len(self))[idx]
        rhs = [list.__getitem__(self, i) for i in indices]
        rv = self.__copy__(shallow=True)
        rv[:] = rhs
        return rv


    def __setitem__(self, idx, a, copy=True):
        """Set idx-th atom to a.

        idx  -- index of atom in this Structure
        a    -- instance of Atom
        copy -- flag for setting to a copy of a.
                When False, set to a and update a.lattice.

        No return value.
        """
        self._uncache('labels')
        adup = copy and Atom(a) or a
        adup.lattice = self.lattice
        list.__setitem__(self, idx, adup)
        return


    def __getslice__(self, lo, hi):
        rv = self.__copy__(shallow=True)
        del rv[hi:], rv[:lo]
        return rv


    def __setslice__(self, lo, hi, atoms, copy=True):
        """Set Structure slice from lo to hi-1 to the sequence of atoms.

        lo    -- low index for the slice
        hi    -- high index of the slice
        atoms -- sequence of Atom instances
        copy  -- flag for using copies of Atom instances.  When False, set
                 to existing instances and update their lattice attributes.

        No return value.
        """
        self._uncache('labels')
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

    # linked atom attributes

    element = _linkAtomAttribute('element',
        '''Character array of atom types.  Assignment updates
        the element attribute of the respective atoms.''',
        toarray=numpy.char.array)

    xyz = _linkAtomAttribute('xyz',
        '''Array of fractional coordinates of all atoms.
        Assignment updates xyz attribute of all atoms.''')

    x = _linkAtomAttribute('x',
        '''Array of all fractional coordinates x.
        Assignment updates xyz attribute of all atoms.''')

    y = _linkAtomAttribute('y',
        '''Array of all fractional coordinates y.
        Assignment updates xyz attribute of all atoms.''')

    z = _linkAtomAttribute('z',
        '''Array of all fractional coordinates z.
        Assignment updates xyz attribute of all atoms.''')

    names = _linkAtomAttribute('name',
        '''Character array of atom names.  Assignment updates
        the name attribute of all atoms.''',
        toarray=numpy.char.array)

    occupancy = _linkAtomAttribute('occupancy',
        '''Array of atom occupancies.  Assignment updates the
        occupancy attribute of all atoms.''')

    xyz_cartn = _linkAtomAttribute('xyz_cartn',
        '''Array of absolute Cartesian coordinates of all atoms.
        Assignment updates the xyz attribute of all atoms.''')

    anisotropy = _linkAtomAttribute('anisotropy',
        '''Boolean array for anisotropic thermal displacement flags.
        Assignment updates the anisotropy attribute of all atoms.''')

    U = _linkAtomAttribute('U',
        '''Array of anisotropic thermal displacement tensors.
        Assignment updates the U and anisotropy attributes of all atoms.''')

    Uisoequiv = _linkAtomAttribute('Uisoequiv',
        '''Array of isotropic thermal displacement or equivalent values.
        Assignment updates the U attribute of all atoms.''')

    U11 = _linkAtomAttribute('U11',
        '''Array of U11 elements of the anisotropic displacement tensors.
        Assignment updates the U and anisotropy attributes of all atoms.''')

    U22 = _linkAtomAttribute('U22',
        '''Array of U22 elements of the anisotropic displacement tensors.
        Assignment updates the U and anisotropy attributes of all atoms.''')

    U33 = _linkAtomAttribute('U33',
        '''Array of U33 elements of the anisotropic displacement tensors.
        Assignment updates the U and anisotropy attributes of all atoms.''')

    U12 = _linkAtomAttribute('U12',
        '''Array of U12 elements of the anisotropic displacement tensors.
        Assignment updates the U and anisotropy attributes of all atoms.''')

    U13 = _linkAtomAttribute('U13',
        '''Array of U13 elements of the anisotropic displacement tensors.
        Assignment updates the U and anisotropy attributes of all atoms.''')

    U23 = _linkAtomAttribute('U23',
        '''Array of U23 elements of the anisotropic displacement tensors.
        Assignment updates the U and anisotropy attributes of all atoms.''')

    Bisoequiv = _linkAtomAttribute('Bisoequiv',
        '''Array of Debye-Waller isotropic thermal displacement or equivalent
        values.  Assignment updates the U attribute of all atoms.''')

    B11 = _linkAtomAttribute('B11',
        '''Array of B11 elements of the Debye-Waller displacement tensors.
        Assignment updates the U and anisotropy attributes of all atoms.''')

    B22 = _linkAtomAttribute('B22',
        '''Array of B22 elements of the Debye-Waller displacement tensors.
        Assignment updates the U and anisotropy attributes of all atoms.''')

    B33 = _linkAtomAttribute('B33',
        '''Array of B33 elements of the Debye-Waller displacement tensors.
        Assignment updates the U and anisotropy attributes of all atoms.''')

    B12 = _linkAtomAttribute('B12',
        '''Array of B12 elements of the Debye-Waller displacement tensors.
        Assignment updates the U and anisotropy attributes of all atoms.''')

    B13 = _linkAtomAttribute('B13',
        '''Array of B13 elements of the Debye-Waller displacement tensors.
        Assignment updates the U and anisotropy attributes of all atoms.''')

    B23 = _linkAtomAttribute('B23',
        '''Array of B23 elements of the Debye-Waller displacement tensors.
        Assignment updates the U and anisotropy attributes of all atoms.''')

    ####################################################################
    # protected methods
    ####################################################################

    def _update_labels(self):
        """Update the _labels dictionary of unique string labels of atoms.

        No return value.
        """
        kv = zip(self.getLabels(), self[:])
        self._labels = dict(kv)
        self._labels_cached = True
        return


    def _uncache(self, *args):
        """Reset cached flag for a list of internal attributes.

        *args -- list of strings, currently supported are "labels"

        No return value.
        Raise AttributeError for any invalid args.
        """
        for a in args:
            attrname = "_" + a + "_cached"
            setattr(self, attrname, False)
        return


# End of class Structure

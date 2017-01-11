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

"""This module defines class Structure.
"""

import collections
import copy
import numpy
from diffpy.Structure.lattice import Lattice
from diffpy.Structure.atom import Atom
from diffpy.Structure.utils import _linkAtomAttribute, atomBareSymbol

##############################################################################
class Structure(list):
    """Structure --> group of atoms

    Structure class is inherited from Python list.  It contains
    a list of Atom instances.  Structure overloads setitem and setslice
    methods so that the lattice attribute of atoms get set to lattice.

    Data members:
        title   -- structure description
        lattice -- coordinate system (instance of Lattice)
        pdffit  -- None or a dictionary of PDFFit-related metadata
    """

    # default values for instance attributes

    title = ''
    _lattice = None
    pdffit = None

    def __init__(self, atoms=None, lattice=None, title=None,
            filename=None, format=None):
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
        # if filename is specified load it and return
        if filename is not None:
            if any((atoms, lattice, title)):
                emsg = "Cannot use filename and atoms arguments together."
                raise ValueError(emsg)
            readkwargs = (format is not None) and {'format' : format} or {}
            self.read(filename, **readkwargs)
            return
        # copy initialization, must be first to allow lattice, title override
        if isinstance(atoms, Structure):
            Structure.__copy__(atoms, self)
        # assign arguments:
        if title is not None:
            self.title = title
        if lattice is not None:
            self.lattice = lattice
        elif self.lattice is None:
            self.lattice = Lattice()
        # insert atoms unless already done by __copy__
        if not len(self) and atoms is not None:
            self.extend(atoms)
        return


    def copy(self):
        '''Return a deep copy of this Structure object.
        '''
        return copy.copy(self)


    def __copy__(self, target=None):
        '''Create a deep copy of this instance.

        target   -- optional target instance for copying, useful for
                    copying a derived class.  Defaults to new instance
                    of the same type as self.

        Return a duplicate instance of this object.
        '''
        if target is None:
            target = Structure()
        elif target is self:
            return target
        # copy attributes as appropriate:
        target.title = self.title
        target.lattice = Lattice(self.lattice)
        target.pdffit = copy.deepcopy(self.pdffit)
        # copy all atoms to the target
        target[:] = self
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
        self.append(a, copy=False)
        return


    def getLastAtom(self):
        """Return Reference to the last Atom in this structure.
        """
        last_atom = self[-1]
        return last_atom


    def assignUniqueLabels(self):
        """Set a unique label string for each atom in this structure.
        The label strings are formatted as "%(baresymbol)s%(index)i",
        where baresymbol is the element right-stripped of "[0-9][+-]".

        No return value.
        """
        elnum = {}
        # support duplicate atom instances
        islabeled = set()
        for a in self:
            if a in islabeled:  continue
            baresmbl = atomBareSymbol(a.element)
            elnum[baresmbl] = elnum.get(baresmbl, 0) + 1
            a.label = baresmbl + str(elnum[baresmbl])
            islabeled.add(a)
        return


    def distance(self, aid0, aid1):
        """Distance between 2 atoms, no periodic boundary conditions.

        aid0 -- zero based index of the first atom or a string label
                such as "Na1"
        aid1 -- zero based index or string label of the second atom.

        Return float.
        Raise IndexError for invalid arguments.
        """
        # lookup by labels
        a0, a1 = self[aid0, aid1]
        return self.lattice.dist(a0.xyz, a1.xyz)


    def angle(self, aid0, aid1, aid2):
        """The bond angle at the second of three atoms in degrees.

        aid0 -- zero based index of the first atom or a string label
                such as "Na1"
        aid1 -- index or string label for the second atom,
                where the angle is formed
        aid2 -- index or string label for the third atom

        Return float.
        Raise IndexError for invalid arguments.
        """
        a0, a1, a2 = self[aid0, aid1, aid2]
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


    def tolist(self):
        '''Return atoms in this Structure as a standard Python list.
        '''
        rv = [a for a in self]
        return rv

    # Overloaded list Methods and Operators ----------------------------------

    def append(self, a, copy=True):
        """Append atom to a structure and update its lattice attribute.

        a    -- instance of Atom
        copy -- flag for appending a copy of a.
                When False, append a and update a.lattice.

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
        if copy:    adups = [Atom(a) for a in atoms]
        else:       adups = atoms
        for a in adups: a.lattice = self.lattice
        list.extend(self, adups)
        return


    def __getitem__(self, idx):
        """Get one or more atoms in this structure.

        idx  -- atom identifier, which can be integer, string or iterable.
                When integer use standard list lookup.  For iterables use
                numpy lookup, this supports integer or boolean flag arrays.
                For string or string-containing iterables lookup the atoms
                by string label.

        Return an Atom instance for integer or string index or a substructure
        in all other cases.  Raise IndexError for invalid index or for
        non-unique atom label.

        Examples:

        stru[0]  -->  first atom in the Structure
        stru[stru.element == 'Na']  -->  substructure of all Na atoms
        stru['Na3']  -->  atom with a unique label 'Na3'
        stru['Na3', 2, 'Cl2']  -->  substructure of three atoms, lookup by
            label is more efficient when done for several atoms at once.
        """
        try:
            value = list.__getitem__(self, idx)
            rv = value
            if type(idx) is slice:
                rv = self.__emptySharedStructure()
                rv.extend(value, copy=False)
            return rv
        except TypeError:
            pass
        # check if there is any string label that should be resolved
        scalarstringlabel = isinstance(idx, basestring)
        hasstringlabel = scalarstringlabel or (
            isinstance(idx, collections.Iterable) and
            any(isinstance(ii, basestring) for ii in idx))
        # if not, use numpy indexing to resolve idx
        if not hasstringlabel:
            idx1 = idx
            if type(idx) is tuple:
                idx1 = numpy.r_[idx]
            indices = numpy.arange(len(self))[idx1]
            rhs = [list.__getitem__(self, i) for i in indices]
            rv = self.__emptySharedStructure()
            rv.extend(rhs, copy=False)
            return rv
        # here we need to resolve at least one string label
        # build a map of labels to indices and mark duplicate labels
        duplicate = object()
        labeltoindex = {}
        for i, a in enumerate(self):
            labeltoindex[a.label] = (
                duplicate if a.label in labeltoindex else i)
        def _resolveindex(aid):
            aid1 = aid
            if type(aid) is str:
                aid1 = labeltoindex.get(aid, None)
                if aid1 is None:
                    raise IndexError("Invalid atom label %r." % aid)
                if aid1 is duplicate:
                    raise IndexError("Atom label %r is not unique." % aid)
            return aid1
        # generate new index object that has no strings
        if scalarstringlabel:
            idx2 = _resolveindex(idx)
        # for iterables preserved the tuple object type
        else:
            idx2 = map(_resolveindex, idx)
            if type(idx) is tuple:
                idx2 = tuple(idx2)
        # call this function again and hope there is no recursion loop
        rv = self[idx2]
        return rv


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


    def __getslice__(self, lo, hi):
        '''Get a slice of atoms from this Structure.

        lo, hi   -- slice indices, negative values are not supported

        Return a sub-structure with atom instances in the slice.
        '''
        rv = self.__emptySharedStructure()
        rv.extend(list.__getslice__(self, lo, hi), copy=False)
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
        if copy:
            ownatoms = set(list.__getslice__(self, lo, hi))
            adups = [(a in ownatoms and a or Atom(a)) for a in atoms]
        else:
            adups = atoms
        for a in adups: a.lattice = self.lattice
        list.__setslice__(self, lo, hi, adups)
        return


    def __add__(self, other):
        '''Return new Structure object with appended atoms from other.

        other    -- sequence of Atom instances

        Return new Structure with a copy of Atom instances.
        '''
        rv = copy.copy(self)
        rv += other
        return rv


    def __iadd__(self, other):
        '''Extend this Structure with atoms from other.

        other    -- sequence of Atom instances

        Return self.
        '''
        self.extend(other)
        return self


    def __sub__(self, other):
        '''Return new Structure that has atoms from the other removed.

        other    -- sequence of Atom instances

        Return new Structure with a copy of Atom instances.
        '''
        otherset = set(other)
        keepindices = [i for i, a in enumerate(self) if not a in otherset]
        rv = copy.copy(self[keepindices])
        return rv


    def __isub__(self, other):
        '''Remove other atoms if present in this structure.

        other    -- sequence of Atom instances

        Return self.
        '''
        otherset = set(other)
        self[:] = [a for a in self if a not in otherset]
        return self


    def __mul__(self, n):
        '''Return new Structure with n-times concatenated atoms from self.
        Atoms and lattice in the new structure are all copies.

        n    -- integer multiple

        Return new Structure.
        '''
        rv = copy.copy(self[:0])
        rv += n * self.tolist()
        return rv

    # right-side multiplication is the same as left-side
    __rmul__ = __mul__


    def __imul__(self, n):
        '''Concatenate this Structure to n-times more atoms.
        For positive multiple the current Atom objects remain at the
        beginning of this Structure.

        n    -- integer multiple

        Return self.
        '''
        if n <= 0:
            self[:] = []
        else:
            self.extend((n - 1) * self.tolist())
        return self

    # Properties -------------------------------------------------------------

    # lattice

    def _get_lattice(self):
        return self._lattice

    def _set_lattice(self, value):
        for a in self:  a.lattice = value
        self._lattice = value
        return

    lattice = property(_get_lattice, _set_lattice, doc =
        "Coordinate system for this Structure.")

    # composition

    def _get_composition(self):
        rv = {}
        for a in self:
            rv[a.element] = rv.get(a.element, 0.0) + a.occupancy
        return rv

    composition = property(_get_composition,
        doc="Dictionary of chemical symbols and their total occupancies.")

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

    label = _linkAtomAttribute('label',
        '''Character array of atom names.  Assignment updates
        the label attribute of all atoms.''',
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

    # Private Methods --------------------------------------------------------

    def __emptySharedStructure(self):
        '''Return empty Structure with standard attributes same as in self.
        '''
        rv = Structure()
        rv.__dict__.update([(k, getattr(self, k)) for k in rv.__dict__])
        return rv

# End of class Structure

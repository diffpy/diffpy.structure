#!/usr/bin/env python
##############################################################################
#
# diffpy.structure  by DANSE Diffraction group
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
"""This module defines class `Structure`."""

import copy as copymod

import numpy
from ase import Atoms as ASEAtoms

from diffpy.structure.atom import Atom
from diffpy.structure.lattice import Lattice
from diffpy.structure.utils import _linkAtomAttribute, atomBareSymbol, isiterable
from diffpy.utils._deprecator import build_deprecation_message, deprecated

# ----------------------------------------------------------------------------

base = "diffpy.structure.Structure"
removal_version = "4.0.0"
assignUniqueLabels_deprecation_msg = build_deprecation_message(
    base,
    "assignUniqueLabels",
    "assign_unique_labels",
    removal_version,
)


class Structure(list):
    """Define group of atoms in a specified lattice. Structure --> group
    of atoms.

    `Structure` class is inherited from Python `list`. It contains
    a list of `Atom` instances. `Structure` overloads `setitem` and `setslice`
    methods so that the `lattice` attribute of atoms get set to `lattice`.

    Parameters
    ----------
    atoms : list of Atom or Structure, Optional
        List of `Atom` instances to be included in this `Structure`.
        When `atoms` argument is an existing `Structure` instance,
        the new structure is its copy.
    lattice : Lattice, Optional
        Instance of `Lattice` defining coordinate systems, property.
    title : str, Optional
        String description of the structure.
    filename : str, Optional
        Name of a file to load the structure from.
    format : str, Optional
        `Structure` format of the loaded `filename`. By default
        all structure formats are tried one by one. Ignored when
        `filename` has not been specified.

    Note
    ----
    Cannot use `filename` and `atoms` arguments together. Overrides `atoms` argument
    when `filename` is specified.

    Attributes
    ----------
    title : str
        String description of the structure.
    lattice : Lattice
        Instance of `Lattice` defining coordinate systems.
    pdffit : None or dict
        Dictionary of PDFFit-related metadata.

    Examples
    --------
    ``Structure(stru)`` create a copy of `Structure` instance stru.

    >>> stru = Structure()
    >>> copystru = Structure(stru)

    `Structure` is inherited from a list it can use list expansions.

    >>> oxygen_atoms = [a for a in stru if a.element == "O" ]
    >>> oxygen_stru = Structure(oxygen_atoms, lattice=stru.lattice)
    """

    # default values for instance attributes
    title = ""
    """str: default values for `title`."""

    _lattice = None
    pdffit = None
    """None: default values for `pdffit`."""

    def __init__(self, atoms=None, lattice=None, title=None, filename=None, format=None):
        # if filename is specified load it and return
        if filename is not None:
            if any((atoms, lattice, title)):
                emsg = "Cannot use filename and atoms arguments together."
                raise ValueError(emsg)
            readkwargs = (format is not None) and {"format": format} or {}
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
        """Return a copy of this `Structure` object."""
        return copymod.copy(self)

    def __copy__(self, target=None):
        """Create a deep copy of this instance.

        Parameters
        ----------
        target :
            Optional target instance for copying, useful for
            copying a derived class. Defaults to new instance
            of the same type as self.

        Returns
        -------
        A duplicate instance of this object.
        """
        if target is None:
            target = Structure()
        elif target is self:
            return target
        # copy attributes as appropriate:
        target.title = self.title
        target.lattice = Lattice(self.lattice)
        target.pdffit = copymod.deepcopy(self.pdffit)
        # copy all atoms to the target
        target[:] = self
        return target

    def __str__(self):
        """Simple string representation."""
        s_lattice = "lattice=%s" % self.lattice
        s_atoms = "\n".join([str(a) for a in self])
        return s_lattice + "\n" + s_atoms

    def addNewAtom(self, *args, **kwargs):
        """Add new `Atom` instance to the end of this `Structure`.

        Parameters
        ----------
        *args, **kwargs :
            See `Atom` class constructor.
        """
        kwargs["lattice"] = self.lattice
        a = Atom(*args, **kwargs)
        self.append(a, copy=False)
        return

    def getLastAtom(self):
        """Return Reference to the last `Atom` in this structure."""
        last_atom = self[-1]
        return last_atom

    def get_chemical_symbols(self, include_charge_state=False):
        """Return list of chemical symbols for all `Atoms` in this
        structure.

        Parameters
        ----------
        include_charge_state : bool, optional
            If ``True``, include charge state in the chemical symbol (e.g., "Fe2+").
        Returns
        -------
        list of str
            The list of chemical symbols for all `Atoms` in this structure.
        """
        symbols_with_charge = [a.element for a in self]
        if include_charge_state:
            return symbols_with_charge
        else:
            symbols = [atomBareSymbol(sym) for sym in symbols_with_charge]
            return symbols

    def get_fractional_coordinates(self):
        """Return array of fractional coordinates of all `Atoms` in this
        structure.

        Returns
        -------
        numpy.ndarray
            The array of fractional coordinates of all `Atoms` in this structure.
        """
        coords = numpy.array([a.xyz for a in self])
        return coords

    def get_cartesian_coordinates(self):
        """Return array of Cartesian coordinates of all `Atoms` in this
        structure.

        Returns
        -------
        numpy.ndarray
            The array of Cartesian coordinates of all `Atoms` in this structure.
        """
        cartn_coords = numpy.array([a.xyz_cartn for a in self])
        return cartn_coords

    def get_anisotropic_displacement_parameters(self):
        """Return array of anisotropic displacement parameters of all
        `Atoms` in this structure.

        Returns
        -------
        numpy.ndarray
            The array of anisotropic displacement parameters of all `Atoms` in this structure.
        """
        adps = numpy.array([a.U for a in self])
        return adps

    def get_isotropic_displacement_parameters(self):
        """Return array of isotropic displacement parameters of all
        `Atoms` in this structure.

        Returns
        -------
        numpy.ndarray
            The array of isotropic displacement parameters of all `Atoms` in this structure.
        """
        idps = numpy.array([a.Uisoequiv for a in self])
        return idps

    def get_occupancies(self):
        """Return array of occupancies of all `Atoms` in this structure.

        Returns
        -------
        numpy.ndarray
            The array of occupancies of all `Atoms` in this structure.
        """
        occupancies = numpy.array([a.occupancy for a in self])
        return occupancies

    def convert_ase_to_diffpy_structure(
        self,
        ase_atoms: ASEAtoms,
        lost_info: list[str] | None = None,
    ) -> Structure | tuple[Structure, dict]:  # noqa
        """Convert ASE `Atoms` object to this `Structure` instance.

        Parameters
        ----------
        ase_structure : ase.Atoms
            The ASE `Atoms` object to be converted.
        lost_info : list of str, optional
            The list of attribute names to extract from the ASE `Atoms`
            object that do not have a direct equivalent in the `Structure` class.
            object that is not currently available in the `Structure` class.
             Default is False.

        Returns
        -------
        Structure
            Reference to this `Structure` object with updated attributes and `Atom` instances.
        lost_info : dict, optional
            The dictionary containing any information from the ASE `Atoms`
            object that is not currently available in the `Structure` class.
            Default behavior is to return only the `Structure` instance.
            If `lost_info` is provided, it will be a dictionary containing
            any information from the ASE `Atoms`.
            This may include information such as magnetic moments, charge states,
            or other ASE-specific properties that do not have a direct equivalent
            in the `Structure` class.

        Raises
        ------
        TypeError
            If the input `ase_structure` is not an instance of `ase.Atoms`.
        ValueError
            If any of the specified `lost_info` attributes are not present in the ASE `Atoms` object.

        Examples
        --------
        An example of converting an `ASE.Atoms` instance to a `Structure` instance,

        .. code-block:: python
            from ase import Atoms
            from diffpy.structure import Structure

            # Create an ASE Atoms object
            ase_atoms = Atoms('H2O', positions=[[0, 0, 0], [0, 0, 1], [1, 0, 0]])

            # Convert to a diffpy Structure object
            structure = Structure()
            structure.convert_ase_to_diffpy(ase_atoms,


        To extract additional information from the ASE `Atoms` object that is not
        directly represented in the `Structure` class, such as magnetic moments,
        you can specify an attribute or method of `ASE.Atoms` as
        a list of strings in `lost_info` list. For example,

        .. code-block:: python
                lost_info = structure.convert_ase_to_diffpy(
                                                            ase_atoms,
                                                            lost_info=['get_magnetic_moments']
                                                            )

        will return a dictionary with the magnetic moments of the atoms in the ASE `Atoms` object.
        """
        if not isinstance(ase_atoms, ASEAtoms):
            raise TypeError("Input must be an instance of ase.Atoms.")
        # --- structure conversion ---
        symbols = ase_atoms.get_chemical_symbols()
        scaled_positions = ase_atoms.get_scaled_positions()
        for sym, xyz in zip(symbols, scaled_positions):
            self.append(Atom(sym, xyz=xyz))
        # --- optional extraction ---
        if lost_info is None:
            return
        extracted_info = {}
        for name in lost_info:
            if not hasattr(ase_atoms, name):
                raise ValueError(f"ASE.Atoms object has no attribute '{name}'.")
            try:
                attr = getattr(ase_atoms, name)
                value = attr() if callable(attr) else attr
                # try to copy (safe for numpy arrays, dicts, etc.)
                try:
                    value = value.copy()
                except Exception:
                    pass
                extracted_info[name] = value
            except Exception as e:
                extracted_info[name] = f"ERROR: {type(e).__name__}: {e}"
        return extracted_info

    def assign_unique_labels(self):
        """Set a unique label string for each `Atom` in this structure.

        The label strings are formatted as "%(baresymbol)s%(index)i",
        where baresymbol is the element right-stripped of "[0-9][+-]".
        """
        elnum = {}
        # support duplicate atom instances
        islabeled = set()
        for a in self:
            if a in islabeled:
                continue
            baresmbl = atomBareSymbol(a.element)
            elnum[baresmbl] = elnum.get(baresmbl, 0) + 1
            a.label = baresmbl + str(elnum[baresmbl])
            islabeled.add(a)
        return

    @deprecated(assignUniqueLabels_deprecation_msg)
    def assignUniqueLabels(self):
        """This function has been deprecated and will be removed in
        version 4.0.0.

        Please use diffpy.structure.Structure.assign_unique_labels
        instead.
        """
        return self.assign_unique_labels()

    def distance(self, aid0, aid1):
        """Calculate distance between 2 `Atoms`, no periodic boundary
        conditions.

        Parameters
        ----------
        aid0 : int or str
            Zero based index of the first `Atom` or a string label.
        aid1 : int or str
            Zero based index or string label of the second atom.

        Returns
        -------
        float
            Distance between the two `Atoms` in Angstroms.

        Raises
        ------
        IndexError
            If any of the `Atom` indices or labels are invalid.
        """
        # lookup by labels
        a0, a1 = self[aid0, aid1]
        return self.lattice.dist(a0.xyz, a1.xyz)

    def angle(self, aid0, aid1, aid2):
        """The bond angle at the second of three `Atoms` in degrees.

        Parameters
        ----------
        aid0 : int or str
            Zero based index of the first `Atom` or a string label.
        aid1 : int or str
            Index or string label for the second atom, where the angle is formed.
        aid2 : int or str
            Index or string label for the third atom.

        Returns
        -------
        float
            The bond angle in degrees.

        Raises
        ------
        IndexError
            If any of the arguments are invalid.
        """
        a0, a1, a2 = self[aid0, aid1, aid2]
        u10 = a0.xyz - a1.xyz
        u12 = a2.xyz - a1.xyz
        return self.lattice.angle(u10, u12)

    def placeInLattice(self, new_lattice):
        """Place structure into `new_lattice` coordinate system.

        Sets `lattice` to `new_lattice` and recalculate fractional coordinates
        of all `Atoms` so their absolute positions remain the same.

        Parameters
        ----------
        new_lattice : Lattice
            New `lattice` to place the structure into.

        Returns
        -------
        Structure
            Reference to this `Structure` object. The `lattice` attribute
            is updated to `new_lattice`.
        """
        Tx = numpy.dot(self.lattice.base, new_lattice.recbase)
        Tu = numpy.dot(self.lattice.normbase, new_lattice.recnormbase)
        for a in self:
            a.xyz = numpy.dot(a.xyz, Tx)
            if a.anisotropy:
                a.U = numpy.dot(numpy.transpose(Tu), numpy.dot(a.U, Tu))
        self.lattice = new_lattice
        return self

    def read(self, filename, format="auto"):
        """Load structure from a file, any original data become lost.

        Parameters
        ----------
        filename : str
            File to be loaded.
        format : str, Optional
            All structure formats are defined in parsers submodule,
            when ``format == 'auto'`` all parsers are tried one by one.

        Returns
        -------
        Parser
            Return instance of data Parser used to process input string. This
            can be inspected for information related to particular format.
        """
        import diffpy.structure
        import diffpy.structure.parsers

        getParser = diffpy.structure.parsers.getParser
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

    def readStr(self, s, format="auto"):
        """Read structure from a string.

        Parameters
        ----------
        s : str
            String with structure definition.
        format : str, Optional
            All structure formats are defined in parsers submodule. When ``format == 'auto'``,
            all parsers are tried one by one.

        Returns
        -------
        Parser
            Return instance of data Parser used to process input string. This
            can be inspected for information related to particular format.
        """
        from diffpy.structure.parsers import getParser

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
        """Save structure to file in the specified format.

        Parameters
        ----------
        filename : str
            File to save the structure to.
        format : str
            `Structure` format to use for saving.

        Note
        ----
        Available structure formats can be obtained by:

            ``from parsers import formats``
        """
        from diffpy.structure.parsers import getParser

        p = getParser(format)
        p.filename = filename
        s = p.tostring(self)
        with open(filename, "w", encoding="utf-8", newline="") as fp:
            fp.write(s)
        return

    def writeStr(self, format):
        """Return string representation of the structure in specified
        format.

        Note
        ----
        Available structure formats can be obtained by:

            ``from parsers import formats``
        """
        from diffpy.structure.parsers import getParser

        p = getParser(format)
        s = p.tostring(self)
        return s

    def tolist(self):
        """Return `Atoms` in this `Structure` as a standard Python
        list."""
        rv = [a for a in self]
        return rv

    # Overloaded list Methods and Operators ----------------------------------

    def append(self, a, copy=True):
        """Append `Atom` to a structure and update its `lattice`
        attribute.

        Parameters
        ----------
        a : Atom
            Instance of `Atom` to be appended.
        copy : bool, Optional
            Flag for appending a copy of `a`. When ``False``, append `a` and update `a.lattice`.
        """
        adup = copy and Atom(a) or a
        adup.lattice = self.lattice
        super(Structure, self).append(adup)
        return

    def insert(self, idx, a, copy=True):
        """Insert `Atom` a before position idx in this `Structure`.

        Parameters
        ----------
        idx : int
            Position in `Atom` list.
        a : Atom
            Instance of `Atom` to be inserted.
        copy : bool, Optional
            Flag for inserting a copy of `a`. When ``False``, append `a` and update `a.lattice`.
        """
        adup = copy and copymod.copy(a) or a
        adup.lattice = self.lattice
        super(Structure, self).insert(idx, adup)
        return

    def extend(self, atoms, copy=None):
        """Extend `Structure` with an iterable of `atoms`.

        Update the `lattice` attribute of all added `atoms`.

        Parameters
        ----------
        atoms : Iterable
            The `Atom` objects to be appended to this `Structure`.
        copy : bool, Optional
            Flag for adding copies of `Atom` objects.
            Make copies when ``True``, append `atoms` unchanged when ``False``.
            The default behavior is to make copies when `atoms` are of
            `Structure` type or if new atoms introduce repeated objects.
        """
        adups = (copymod.copy(a) for a in atoms)
        if copy is None:
            if isinstance(atoms, Structure):
                newatoms = adups
            else:
                memo = set(id(a) for a in self)

                def nextatom(a):
                    return a if id(a) not in memo else copymod.copy(a)

                def mark(a):
                    return (memo.add(id(a)), a)[-1]

                newatoms = (mark(nextatom(a)) for a in atoms)
        elif copy:
            newatoms = adups
        else:
            newatoms = atoms

        def setlat(a):
            return (setattr(a, "lattice", self.lattice), a)[-1]

        super(Structure, self).extend(setlat(a) for a in newatoms)
        return

    def __getitem__(self, idx):
        """Get one or more `Atoms` in this structure.

        Parameters
        ----------
        idx : int or str or Iterable
            `Atom` identifier. When integer use standard list lookup.
            For iterables use numpy lookup, this supports integer or
            boolean flag arrays. For string or string-containing iterables
            lookup the `Atoms` by string label.

        Returns
        -------
        Atom or Structure
            An `Atom` instance for integer or string index or a substructure
            in all other cases.

        Raises
        ------
        IndexError
            If the index is invalid or the `Atom` label is not unique.

        Examples
        --------
        First `Atom` in the `Structure`:

        >>> stru[0]

        Substructure of all ``'Na'`` `Atoms`:

        >>> stru[stru.element == 'Na']

        `Atom` with a unique label ``'Na3'``:
        >>> stru['Na3']

        Substructure of three `Atoms`, lookup by label is more efficient
        when done for several `Atoms` at once.

        >>> stru['Na3', 2, 'Cl2']
        """
        if isinstance(idx, slice):
            rv = self.__emptySharedStructure()
            lst = super(Structure, self).__getitem__(idx)
            rv.extend(lst, copy=False)
            return rv
        try:
            rv = super(Structure, self).__getitem__(idx)
            return rv
        except TypeError:
            pass
        # check if there is any string label that should be resolved
        scalarstringlabel = isinstance(idx, str)
        hasstringlabel = scalarstringlabel or (isiterable(idx) and any(isinstance(ii, str) for ii in idx))
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
            labeltoindex[a.label] = duplicate if a.label in labeltoindex else i

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
        # for iterables preserve the tuple object type
        else:
            idx2 = [_resolveindex(i) for i in idx]
            if type(idx) is tuple:
                idx2 = tuple(idx2)
        # call this function again and hope there is no recursion loop
        rv = self[idx2]
        return rv

    def __setitem__(self, idx, value, copy=True):
        """Assign `self[idx]` `Atom` to value.

        Parameters
        ----------
        idx : int or slice
            Index of `Atom` in this `Structure` or a slice.
        value : Atom or Iterable
            Instance of `Atom` or an iterable.
        copy : bool, Optional
            Flag for making a copy of the value. When ``False``, update
            the `lattice` attribute of `Atom` objects present in value.
            Default is ``True``.
        """
        # handle slice assignment
        if isinstance(idx, slice):

            def _fixlat(a):
                a.lattice = self.lattice
                return a

            v1 = value
            if copy:
                keep = set(super(Structure, self).__getitem__(idx))
                v1 = (a if a in keep else Atom(a) for a in value)
            vfinal = filter(_fixlat, v1)
        # handle scalar assignment
        else:
            vfinal = Atom(value) if copy else value
            vfinal.lattice = self.lattice
        super(Structure, self).__setitem__(idx, vfinal)
        return

    def __add__(self, other):
        """Return new `Structure` object with appended `Atoms` from
        other.

        Parameters
        ----------
        other : sequence of Atom
            Sequence of `Atom` instances.

        Returns
        -------
        Structure
            New `Structure` with a copy of `Atom` instances.
        """
        rv = copymod.copy(self)
        rv += other
        return rv

    def __iadd__(self, other):
        """Extend this `Structure` with `Atoms` from other.

        Parameters
        ----------
        other : sequence of Atom
            Sequence of `Atom` instances.

        Returns
        -------
        Structure
            Reference to this `Structure` object.
        """
        self.extend(other, copy=True)
        return self

    def __sub__(self, other):
        """Return new `Structure` that has `Atoms` from the other
        removed.

        Parameters
        ----------
        other : sequence of Atom
            Sequence of `Atom` instances.

        Returns
        -------
        Structure
            New `Structure` with a copy of `Atom` instances.
        """
        otherset = set(other)
        keepindices = [i for i, a in enumerate(self) if a not in otherset]
        rv = copymod.copy(self[keepindices])
        return rv

    def __isub__(self, other):
        """Remove other `Atoms` if present in this structure.

        Parameters
        ----------
        other : sequence of Atom
            Sequence of `Atom` instances.

        Returns
        -------
        Structure
            Reference to this `Structure` object.
        """
        otherset = set(other)
        self[:] = [a for a in self if a not in otherset]
        return self

    def __mul__(self, n):
        """Return new `Structure` with n-times concatenated `Atoms` from
        self. `Atoms` and `lattice` in the new structure are all copies.

        Parameters
        ----------
        n : int
            Integer multiple.

        Returns
        -------
        Structure
            New `Structure` with n-times concatenated `Atoms`.
        """
        rv = copymod.copy(self[:0])
        rv += n * self.tolist()
        return rv

    # right-side multiplication is the same as left-side
    __rmul__ = __mul__

    def __imul__(self, n):
        """Concatenate this `Structure` to n-times more `Atoms`. For
        positive multiple the current `Atom` objects remain at the
        beginning of this `Structure`.

        Parameters
        ----------
        n : int
            Integer multiple.

        Returns
        -------
        Structure
            Reference to this `Structure` object.
        """
        if n <= 0:
            self[:] = []
        else:
            self.extend((n - 1) * self.tolist(), copy=True)
        return self

    # Properties -------------------------------------------------------------

    # lattice

    def _get_lattice(self):
        return self._lattice

    def _set_lattice(self, value):
        for a in self:
            a.lattice = value
        self._lattice = value
        return

    lattice = property(
        _get_lattice,
        _set_lattice,
        doc="Coordinate system for this `Structure`.",
    )

    # composition

    def _get_composition(self):
        rv = {}
        for a in self:
            rv[a.element] = rv.get(a.element, 0.0) + a.occupancy
        return rv

    composition = property(
        _get_composition,
        doc="Dictionary of chemical symbols and their total occupancies.",
    )

    # linked atom attributes

    element = _linkAtomAttribute(
        "element",
        """Character array of `Atom` types. Assignment updates
        the element attribute of the respective `Atoms`.
        Set the maximum length of the element string to 5 characters.""",
        toarray=lambda items: numpy.char.array(items, itemsize=5),
    )

    xyz = _linkAtomAttribute(
        "xyz",
        """Array of fractional coordinates of all `Atoms`.
        Assignment updates `xyz` attribute of all `Atoms`.""",
    )

    x = _linkAtomAttribute(
        "x",
        """Array of all fractional coordinates `x`.
        Assignment updates `xyz` attribute of all `Atoms`.""",
    )

    y = _linkAtomAttribute(
        "y",
        """Array of all fractional coordinates `y`.
        Assignment updates `xyz` attribute of all `Atoms`.""",
    )

    z = _linkAtomAttribute(
        "z",
        """Array of all fractional coordinates `z`.
        Assignment updates `xyz` attribute of all `Atoms`.""",
    )

    label = _linkAtomAttribute(
        "label",
        """Character array of `Atom` names. Assignment updates
        the label attribute of all `Atoms`.
        Set the maximum length of the label string to 5 characters.""",
        toarray=lambda items: numpy.char.array(items, itemsize=5),
    )

    occupancy = _linkAtomAttribute(
        "occupancy",
        """Array of `Atom` occupancies. Assignment updates the
        occupancy attribute of all `Atoms`.""",
    )

    xyz_cartn = _linkAtomAttribute(
        "xyz_cartn",
        """Array of absolute Cartesian coordinates of all `Atoms`.
        Assignment updates the `xyz` attribute of all `Atoms`.""",
    )

    anisotropy = _linkAtomAttribute(
        "anisotropy",
        """Boolean array for anisotropic thermal displacement flags.
        Assignment updates the anisotropy attribute of all `Atoms`.""",
    )

    U = _linkAtomAttribute(
        "U",
        """Array of anisotropic thermal displacement tensors.
        Assignment updates the U and anisotropy attributes of all `Atoms`.""",
    )

    Uisoequiv = _linkAtomAttribute(
        "Uisoequiv",
        """Array of isotropic thermal displacement or equivalent values.
        Assignment updates the U attribute of all `Atoms`.""",
    )

    U11 = _linkAtomAttribute(
        "U11",
        """Array of `U11` elements of the anisotropic displacement tensors.
        Assignment updates the U and anisotropy attributes of all `Atoms`.""",
    )

    U22 = _linkAtomAttribute(
        "U22",
        """Array of `U22` elements of the anisotropic displacement tensors.
        Assignment updates the U and anisotropy attributes of all `Atoms`.""",
    )

    U33 = _linkAtomAttribute(
        "U33",
        """Array of `U33` elements of the anisotropic displacement tensors.
        Assignment updates the U and anisotropy attributes of all `Atoms`.""",
    )

    U12 = _linkAtomAttribute(
        "U12",
        """Array of `U12` elements of the anisotropic displacement tensors.
        Assignment updates the U and anisotropy attributes of all `Atoms`.""",
    )

    U13 = _linkAtomAttribute(
        "U13",
        """Array of `U13` elements of the anisotropic displacement tensors.
        Assignment updates the U and anisotropy attributes of all `Atoms`.""",
    )

    U23 = _linkAtomAttribute(
        "U23",
        """Array of `U23` elements of the anisotropic displacement tensors.
        Assignment updates the U and anisotropy attributes of all `Atoms`.""",
    )

    Bisoequiv = _linkAtomAttribute(
        "Bisoequiv",
        """Array of Debye-Waller isotropic thermal displacement or equivalent
        values. Assignment updates the U attribute of all `Atoms`.""",
    )

    B11 = _linkAtomAttribute(
        "B11",
        """Array of `B11` elements of the Debye-Waller displacement tensors.
        Assignment updates the U and anisotropy attributes of all `Atoms`.""",
    )

    B22 = _linkAtomAttribute(
        "B22",
        """Array of `B22` elements of the Debye-Waller displacement tensors.
        Assignment updates the U and anisotropy attributes of all `Atoms`.""",
    )

    B33 = _linkAtomAttribute(
        "B33",
        """Array of `B33` elements of the Debye-Waller displacement tensors.
        Assignment updates the U and anisotropy attributes of all `Atoms`.""",
    )

    B12 = _linkAtomAttribute(
        "B12",
        """Array of `B12` elements of the Debye-Waller displacement tensors.
        Assignment updates the U and anisotropy attributes of all `Atoms`.""",
    )

    B13 = _linkAtomAttribute(
        "B13",
        """Array of `B13` elements of the Debye-Waller displacement tensors.
        Assignment updates the U and anisotropy attributes of all `Atoms`.""",
    )

    B23 = _linkAtomAttribute(
        "B23",
        """Array of `B23` elements of the Debye-Waller displacement tensors.
        Assignment updates the U and anisotropy attributes of all `Atoms`.""",
    )

    # Private Methods --------------------------------------------------------

    def __emptySharedStructure(self):
        """Return empty `Structure` with standard attributes same as in
        self."""
        rv = Structure()
        rv.__dict__.update([(k, getattr(self, k)) for k in rv.__dict__])
        return rv


# End of class Structure

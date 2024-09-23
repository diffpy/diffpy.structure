#!/usr/bin/env python
##############################################################################
#
# diffpy.structure  by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2006 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################

"""
Provide class Atom for managing properties of an atom in structure model.
"""

import numpy

from diffpy.structure.lattice import cartesian as cartesian_lattice

# conversion constants
_BtoU = 1.0 / (8 * numpy.pi**2)
_UtoB = 1.0 / _BtoU

# ----------------------------------------------------------------------------


class Atom(object):
    """Storage of structure information relevant for a single atom.

    This class manages atom information such as element symbol, position
    in fractional and Cartesian coordinates, atomic displacement parameters
    and so forth.

    Parameters
    ----------
    atype : str or Atom, Optional
        The string atom type to be set as the `element` attribute.
        By default an empty string. When of the `Atom` type, create
        a copy of `atype` and adjust it per other arguments.
    xyz : numpy.ndarray, Optional
        Fractional coordinates within the associated `lattice`.
        By default ``[0, 0, 0]``.
    label : str, Optional
        A unique string `label` for referring to this `Atom`.
        By default an empty string.
    occupancy : float, Optional
        The initial `occupancy` of this atom, by default ``1``.
    anisotropy : bool, Optional
        The flag for anisotropic thermal displacements parameters.
        This overrides `anisotropy` implied by presence of the
        *U* or *Uisoequiv* arguments. Defaults to ``False``
        when not set in any other way.
    U : numpy.ndarray, Optional
        The 3x3 matrix of anisotropic thermal displacement parameters.
        When present `anisotropy` defaults to ``True``.
    Uisoequiv: float, Optional
        The isotropic atomic displacement parameter. The `anisotropy`
        defaults to ``False`` when present. Only one of the *U* and
        *Uisoequiv* arguments may be provided at the same time. Assume
        zero atomic displacements when *U* and *Uisoequiv* are unset.
    lattice : Lattice, Optional
        Coordinate system for the fractional coordinates `xyz`.
        Use the absolute Cartesian system when ``None``.

    Attributes
    ----------
    element : str
        The string type of the atom. An element or ion symbol.
    xyz : numpy.ndarray
        The fractional coordinates in the associated `lattice`.
    label : str
        A unique string label referring to this atom, for example, "C_1".
        The *label* can be used to reference this atom when contained in
        a `Structure` object.
    occupancy : float
        The fractional occupancy of this atom.
    lattice : Lattice
        Coordinate system for the fractional coordinates `xyz` and
        the tensor of atomic displacement parameters `U`.
        Use the absolute Cartesian coordinates when ``None``.

    Note
    ----
    Cannot use both U and Uisoequiv arguments at the same time.
    """

    # Private attributes
    #
    #   _U : 3-by-3 ndarray
    #       Internal storage of the displacement parameters.

    # instance attributes that have immutable default values
    element = ""
    """str: Default values of `element`."""

    label = ""
    """str: Default values of `label`."""

    occupancy = 1.0
    """float: Default values of `occupancy`."""

    _anisotropy = False
    lattice = None
    """None: Default values of `lattice`."""

    def __init__(
        self,
        atype=None,
        xyz=None,
        label=None,
        occupancy=None,
        anisotropy=None,
        U=None,
        Uisoequiv=None,
        lattice=None,
    ):
        # check arguments
        if U is not None and Uisoequiv is not None:
            emsg = "Cannot use both U and Uisoequiv arguments."
            raise ValueError(emsg)
        # declare data members
        self.xyz = numpy.zeros(3, dtype=float)
        self._U = numpy.zeros((3, 3), dtype=float)
        # assign them as needed
        if isinstance(atype, Atom):
            atype.__copy__(target=self)
        elif atype is not None:
            self.element = atype
        # take care of remaining arguments
        if xyz is not None:
            self.xyz[:] = xyz
        if label is not None:
            self.label = label
        if occupancy is not None:
            self.occupancy = float(occupancy)
        if U is not None:
            self.anisotropy = True
            self._U[:] = U
        if Uisoequiv is not None:
            self.anisotropy = False
            self.Uisoequiv = Uisoequiv
        # lattice needs to be set before anisotropy
        if lattice is not None:
            self.lattice = lattice
        # process anisotropy after U, Uisoequiv and lattice.
        if anisotropy is not None:
            self.anisotropy = bool(anisotropy)
        return

    def msdLat(self, vl):
        """Calculate mean square displacement along the lattice vector.

        Parameters
        ----------
        vl : array_like
            The vector in lattice coordinates.

        Returns
        -------
        float
            The mean square displacement along *vl*.
        """
        if not self.anisotropy:
            return self.Uisoequiv
        # here we need to calculate msd
        lat = self.lattice or cartesian_lattice
        vln = numpy.array(vl, dtype=float) / lat.norm(vl)
        G = lat.metrics
        rhs = numpy.array([G[0] * lat.ar, G[1] * lat.br, G[2] * lat.cr], dtype=float)
        rhs = numpy.dot(rhs, vln)
        msd = numpy.dot(rhs, numpy.dot(self.U, rhs))
        return msd

    def msdCart(self, vc):
        """Calculate mean square displacement along the Cartesian vector.

        Parameters
        ----------
        vc : array_like
            Vector in Cartesian coordinates.

        Returns
        -------
        float
            The mean square displacement along *vc*.
        """
        if not self.anisotropy:
            return self.Uisoequiv
        # here we need to calculate msd
        lat = self.lattice or cartesian_lattice
        vcn = numpy.array(vc, dtype=float)
        vcn /= numpy.sqrt(numpy.sum(vcn**2))
        F1 = lat.normbase
        Uc = numpy.dot(numpy.transpose(F1), numpy.dot(self._U, F1))
        msd = numpy.dot(vcn, numpy.dot(Uc, vcn))
        return msd

    def __repr__(self):
        """String representation of this Atom."""
        xyz = self.xyz
        s = "%-4s %8.6f %8.6f %8.6f %6.4f" % (self.element, xyz[0], xyz[1], xyz[2], self.occupancy)
        return s

    def __copy__(self, target=None):
        """Create a copy of this instance.

        Parameters
        ----------
        target : Atom, Optional
            An already existing `Atom` object to be updated to a duplicate
            of this `Atom`. Create a new Atom object when not specified.
            This facilitates extension of the `__copy__` method
            in a derived class.

        Returns
        -------
        Atom
            The copy of this object.
        """
        if target is None:
            target = Atom()
        elif target is self:
            return target
        target.__dict__.update(self.__dict__)
        target.xyz = numpy.copy(self.xyz)
        target._U = numpy.copy(self._U)
        return target

    # property handlers ------------------------------------------------------

    x = property(
        lambda self: self.xyz[0],
        lambda self, val: self.xyz.__setitem__(0, val),
        doc="float : fractional coordinate *x*, same as ``xyz[0]``.",
    )
    y = property(
        lambda self: self.xyz[1],
        lambda self, val: self.xyz.__setitem__(1, val),
        doc="float : fractional coordinate *y*, same as ``xyz[1]``.",
    )
    z = property(
        lambda self: self.xyz[2],
        lambda self, val: self.xyz.__setitem__(2, val),
        doc="float : fractional coordinate *z*, same as ``xyz[2]``.",
    )

    # xyz_cartn

    @property
    def xyz_cartn(self):
        """numpy.ndarray: Atom position in absolute Cartesian coordinates.

        This is computed from fractional coordinates `xyz` and the
        current `lattice` setup. Assignment to *xyz_cartn* or
        its components is applied on fractional coordinates `xyz`.
        """
        if not self.lattice:
            rv = self.xyz
        else:
            rv = _AtomCartesianCoordinates(self)
        return rv

    @xyz_cartn.setter
    def xyz_cartn(self, value):
        if not self.lattice:
            self.xyz[:] = value
        else:
            self.xyz[:] = self.lattice.fractional(value)
        return

    # anisotropy

    @property
    def anisotropy(self):
        """bool : Flag for allowing anisotropic displacement parameters.

        When ``False`` the tensor of thermal displacement parameters `U`
        must be isotropic and only its diagonal elements are taken into
        account.
        """
        return self._anisotropy

    @anisotropy.setter
    def anisotropy(self, value):
        if bool(value) is self._anisotropy:
            return
        # convert from isotropic to anisotropic
        if value:
            self._U = self.U
        # otherwise convert from anisotropic to isotropic
        else:
            self._U[0, 0] = self.Uisoequiv
        self._anisotropy = bool(value)
        return

    # U

    @property
    def U(self):
        """numpy.ndarray : The 3x3 matrix of anisotropic atomic displacements.

        For isotropic displacements (when `anisotropy` is ``False``)
        assignment to *U* uses only the first ``Unew[0, 0]`` element
        and the remaining components of *U* are adjusted to obtain
        isotropic tensor in the active `lattice`.

        Note
        ----
        Elements of the *U* tensor such as ``U[0, 1]`` should be
        considered read-only as setting them directly leads to
        undefined behavior. Use the `U11`, `U22`, ..., or `B11`,
        `B22`, ..., descriptors to set only some *U* components.
        """
        if not self.anisotropy:
            # for isotropic displacements assume first element
            # to be equal to the displacement value
            lat = self.lattice or cartesian_lattice
            numpy.multiply(self._U[0, 0], lat.isotropicunit, out=self._U)
        return self._U

    @U.setter
    def U(self, value):
        self._U[:] = value
        return

    # Uij elements

    def _get_Uij(self, i, j):
        """The getter function for the `U11`, `U22`, ..., properties."""
        if self.anisotropy:
            return self._U[i, j]
        lat = self.lattice or cartesian_lattice
        return self._U[0, 0] * lat.isotropicunit[i, j]

    def _set_Uij(self, i, j, value):
        """The setter function for the `U11`, `U22`, ..., properties."""
        self._U[i, j] = value
        self._U[j, i] = value
        if not self._anisotropy and i == j != 0:
            self._U[0, 0] = value
        return

    # _doc_uii, _doc_uij are temporary local variables.

    _doc_uii = """
        float : The ``U[{0}, {0}]`` component of the displacement tensor `U`.

        When `anisotropy` is ``False`` setting a new value updates entire
        tensor *U*.
        """

    U11 = property(
        lambda self: self._get_Uij(0, 0), lambda self, value: self._set_Uij(0, 0, value), doc=_doc_uii.format(0)
    )
    U22 = property(
        lambda self: self._get_Uij(1, 1), lambda self, value: self._set_Uij(1, 1, value), doc=_doc_uii.format(1)
    )
    U33 = property(
        lambda self: self._get_Uij(2, 2), lambda self, value: self._set_Uij(2, 2, value), doc=_doc_uii.format(2)
    )

    _doc_uij = """
        float : The ``U[{0}, {1}]`` element of the displacement tensor `U`.

        Sets ``U[{1}, {0}]`` together with ``U[{0}, {1}]``. Assignment
        has no effect when `anisotropy` is ``False``.
        """

    U12 = property(
        lambda self: self._get_Uij(0, 1), lambda self, value: self._set_Uij(0, 1, value), doc=_doc_uij.format(0, 1)
    )
    U13 = property(
        lambda self: self._get_Uij(0, 2), lambda self, value: self._set_Uij(0, 2, value), doc=_doc_uij.format(0, 2)
    )
    U23 = property(
        lambda self: self._get_Uij(1, 2), lambda self, value: self._set_Uij(1, 2, value), doc=_doc_uij.format(1, 2)
    )

    # clean local variables
    del _doc_uii, _doc_uij

    # Uisoequiv

    @property
    def Uisoequiv(self):
        """float : The isotropic displacement parameter or an equivalent value.

        Setting a new value rescales tensor `U` so it yields equivalent
        direction-averaged displacements.
        """
        if not self.anisotropy:
            return self._U[0, 0]
        if self.lattice is None:
            return numpy.trace(self._U) / 3.0
        lat = self.lattice
        rv = (
            1.0
            / 3.0
            * (
                self._U[0, 0] * lat.ar * lat.ar * lat.a * lat.a
                + self._U[1, 1] * lat.br * lat.br * lat.b * lat.b
                + self._U[2, 2] * lat.cr * lat.cr * lat.c * lat.c
                + 2 * self._U[0, 1] * lat.ar * lat.br * lat.a * lat.b * lat.cg
                + 2 * self._U[0, 2] * lat.ar * lat.cr * lat.a * lat.c * lat.cb
                + 2 * self._U[1, 2] * lat.br * lat.cr * lat.b * lat.c * lat.ca
            )
        )
        return rv

    @Uisoequiv.setter
    def Uisoequiv(self, value):
        if self.anisotropy:
            lat = self.lattice or cartesian_lattice
            uequiv = self.Uisoequiv
            if abs(uequiv) < lat._epsilon:
                self._U = value * lat.isotropicunit
            else:
                self._U *= value / uequiv
        else:
            self._U[0, 0] = value
        return

    # Bij elements

    # _doc_bii, _doc_bij are local variables.

    _doc_bii = """
        float : The ``B{0}{0}`` element of the Debye-Waller matrix.

        This is equivalent to ``8 * pi**2 * U{0}{0}``. When `anisotropy`
        is ``False`` setting a new value updates entire tensor `U`.
        """

    _doc_bij = """
        float : The ``B{0}{1}`` element of the Debye-Waller matrix.

        This is equivalent to ``8 * pi**2 * U{0}{1}``. Setting a new
        value updates `U` in a symmetric way. Assignment has no effect
        when `anisotropy` is ``False``.
        """

    B11 = property(
        lambda self: _UtoB * self._get_Uij(0, 0),
        lambda self, value: self._set_Uij(0, 0, _BtoU * value),
        doc=_doc_bii.format(1),
    )
    B22 = property(
        lambda self: _UtoB * self._get_Uij(1, 1),
        lambda self, value: self._set_Uij(1, 1, _BtoU * value),
        doc=_doc_bii.format(2),
    )
    B33 = property(
        lambda self: _UtoB * self._get_Uij(2, 2),
        lambda self, value: self._set_Uij(2, 2, _BtoU * value),
        doc=_doc_bii.format(3),
    )
    B12 = property(
        lambda self: _UtoB * self._get_Uij(0, 1),
        lambda self, value: self._set_Uij(0, 1, _BtoU * value),
        doc=_doc_bij.format(1, 2),
    )
    B13 = property(
        lambda self: _UtoB * self._get_Uij(0, 2),
        lambda self, value: self._set_Uij(0, 2, _BtoU * value),
        doc=_doc_bij.format(1, 3),
    )
    B23 = property(
        lambda self: _UtoB * self._get_Uij(1, 2),
        lambda self, value: self._set_Uij(1, 2, _BtoU * value),
        doc=_doc_bij.format(2, 3),
    )

    # clean local variables
    del _doc_bii, _doc_bij

    # Bisoequiv

    @property
    def Bisoequiv(self):
        """float : The Debye-Waller isotropic displacement or an equivalent value.

        This equals ``8 * pi**2 * Uisoequiv``. Setting a new value
        rescales `U` tensor to yield equivalent direction-average of
        Debye-Waller displacements.
        """
        return _UtoB * self.Uisoequiv

    @Bisoequiv.setter
    def Bisoequiv(self, value):
        self.Uisoequiv = _BtoU * value
        return


# End of class Atom

# Local Helpers --------------------------------------------------------------


class _AtomCartesianCoordinates(numpy.ndarray):
    """Specialized `numpy.ndarray` for accessing Cartesian coordinates.

    Inplace assignments to this array are applied on the *xyz* position
    position of owner `Atom` as per the associated `Atom.lattice`.

    Parameters
    ----------
    atom : Atom
        `Atom` instance to be linked to these coordinate array.
    """

    def __new__(self, atom):
        """Create the underlying numpy array base object."""
        return numpy.empty(3, dtype=float).view(self)

    def __init__(self, atom):
        self._atom = atom
        self.asarray[:] = atom.lattice.cartesian(atom.xyz)
        return

    @property
    def asarray(self):
        """ndarray : This array viewed as standard numpy array."""
        return self.view(numpy.ndarray)

    def __setitem__(self, idx, value):
        """Set some element or slice of this Cartesian coordinates.

        This overrides inplace array assignment to update the
        *xyz* fractional coordinate of the linked `Atom`.
        """
        self.asarray[idx] = value
        self._atom.xyz[:] = self._atom.lattice.fractional(self)
        return

    def __array_wrap__(self, out_arr, context=None, return_scalar=None):
        """Ensure math operations on this type yield standard numpy array."""
        return out_arr.view(numpy.ndarray)


# End of _AtomCartesianCoordinates

#!/usr/bin/env python3
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

"""class Atom for storing properties of a single atom"""

import numpy
from diffpy.structure.lattice import cartesian as cartesian_lattice

# conversion constants
_BtoU = 1.0/(8 * numpy.pi**2)
_UtoB = 1.0/_BtoU

# ----------------------------------------------------------------------------

class Atom:
    """
    Storage of structure information relevant for a single atom.

    This class handles atom attributes such as element type, position in
    fractional and Cartesian coordinates, displacement parameters and
    so forth.

    Parameters
    ----------
    atype : str
            Element symbol string or Atom instance.
    xyz : ndarray
        Fractional coordinates within the associated `lattice`.
    label : str
        String atom label.
    occupancy : float
        Fractional occupancy of this atom.
    anisotropy : bool
        flag for anisotropic thermal displacements parameters.
        Overrides anisotropy implied from the presence of
        ``U`` or ``Uisoequiv`` arguments.
    U : ndarray
        Anisotropic thermal displacement tensor.
        Matrix of anisotropic displacement parameters.
        When specified also set ``anisotropy`` = True.
    Uisoequiv: float
        isotropic displacement parameter.
        When specified anisotropy is set to False.
        Only one of the ``U``, ``Uisoequiv`` arguments can be used.
    lattice : Lattice
        Coordinate system for the fractional coordinates `xyz`.
        When ``None`` use the absolute Cartesian system.

    Attributes
    ----------
    element : str
        String type of the atom.  Element or ion symbol.
    xyz : ndarray
        Fractional coordinates within the associated `lattice`.
    label : str
        String label for this atom such as "C_1".  It can be used
        to reference this atom when owned by some `Structure`.
    occupancy : float
        Fractional occupancy of this atom.
    lattice : Lattice
        Coordinate system for the fractional coordinates `xyz`.
        When ``None`` use the absolute Cartesian system.
    """



    # Private attributes
    #
    #   _U : 3-by-3 ndarray
    #       Internal storage of the displacement parameters.

    # instance attributes that have inmutable default values
    element = ''
    label = ''
    occupancy = 1.0
    _anisotropy = False
    lattice = None

    def __init__(self, atype=None, xyz=None, label=None, occupancy=None,
            anisotropy=None, U=None, Uisoequiv=None, lattice=None):
        """
        Create atom of a specified type at given lattice coordinates. Atom(a) creates a copy of Atom instance a.

        """
        # check arguments
        if U is not None and Uisoequiv is not None:
            emsg = "Cannot use both U and Uisoequiv arguments."
            raise ValueError(emsg)
        # declare data members
        self.xyz = numpy.zeros(3, dtype=float)
        self._U = numpy.zeros((3,3), dtype=float)
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
        """
        Mean square displacement of an atom along lattice vector.

        Parameters
        ----------
        vl : ndarray
            Vector in lattice coordinates.

        Returns
        -------
        msd : float
            mean square displacement.
        """
        if not self.anisotropy:     return self.Uisoequiv
        # here we need to calculate msd
        lat = self.lattice or cartesian_lattice
        vln = numpy.array(vl, dtype=float)/lat.norm(vl)
        G = lat.metrics
        rhs = numpy.array([ G[0]*lat.ar,
                            G[1]*lat.br,
                            G[2]*lat.cr ], dtype=float)
        rhs = numpy.dot(rhs, vln)
        msd = numpy.dot(rhs, numpy.dot(self.U, rhs))
        return msd


    def msdCart(self, vc):
        """
        Mean quare displacement of an atom along cartesian vector.

        Parameters
        ----------
        vc : ndarray
            Vector in absolute cartesian coordinates.

        Returns
        -------
        msd : float
            mean square displacement.
        """
        if not self.anisotropy:     return self.Uisoequiv
        # here we need to calculate msd
        lat = self.lattice or cartesian_lattice
        vcn = numpy.array(vc, dtype=float)
        vcn /= numpy.sqrt(numpy.sum(vcn**2))
        F1 = lat.normbase
        Uc = numpy.dot(numpy.transpose(F1), numpy.dot(self._U, F1))
        msd = numpy.dot(vcn, numpy.dot(Uc, vcn))
        return msd


    def __repr__(self):
        """
        String representation of the atom's ``atype``, ``xyz``, and ``occupancy``.

        Parameters
        ----------
        None.

        Returns
        -------
        s : str
            atom's ``atype``, ``xyz``, and ``occupancy``.
        """
        xyz = self.xyz
        s = "%-4s %8.6f %8.6f %8.6f %6.4f" % \
                (self.element, xyz[0], xyz[1], xyz[2], self.occupancy)
        return s


    def __copy__(self, target=None):
        """
        Create a copy of this instance.

        Parameters
        ----------
        target : Atom, optional
            An already existing Atom object to be updated to a duplicate
            of this Atom. Create a new Atom object when not specified.
            This argument facilitates extension of the `__copy__` method
            in a derived class.

        Returns
        -------
        target : Atom
            A copy of atom class from self.
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

    x = property(lambda self: self.xyz[0],
            lambda self, val: self.xyz.__setitem__(0, val),
            doc='float : fractional coordinate x.')
    y = property(lambda self: self.xyz[1],
            lambda self, val: self.xyz.__setitem__(1, val),
            doc='float : fractional coordinate y.')
    z = property(lambda self: self.xyz[2],
            lambda self, val: self.xyz.__setitem__(2, val),
            doc='float : fractional coordinate z.')

    # xyz_cartn

    @property
    def xyz_cartn(self):
        """
        ndarray: Atom position in the absolute Cartesian coordinates.

        This is computed from fractional coordinates `xyz` and the
        current lattice configuration.  Setting *xyz_cartn* or any
        of its items such as ``xyz_cartn[0]`` updates the fractional
        coordinates `xyz`.
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
            self.xyz = self.lattice.fractional(value)
        return

    # anisotropy

    @property
    def anisotropy(self):
        """
        bool: Allow anisotropic thermal displacement parameters ``U`` if True.

        When ``False`` the anisotropic thermal displacement parameters matrix ``U``
        must be isotropic and only its diagonal elements are taken into account.
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

    # TODO: convert to new-style property definition
    @property
    def U(self):
        """
        ndarray: Anisotropic thermal displacement parameters matrix.

        When specified also set ``anisotropy`` = True.
        Obtain the Anisotropic thermal displacement tensor ``U``.
        """
        if not self.anisotropy:
            # for isotropic displacements assume first element
            # to be equal to the displacement value
            lat = self.lattice or cartesian_lattice
            self._U = self._U[0, 0] * lat.isotropicunit
        return self._U

    @U.setter
    def U(self, value):
        self._U[:] = value
        return

    # Uij elements

    def _get_Uij(self, i, j):
        """
        Obtain the elements of anisotropic thermal displacement tensor ``Uij``.

        Parameters
        ----------
        i : int
            The index of anisotropic thermal displacement tensor ``Uij``.
            i is from (1, 2, 3).
        j : int
            The index of anisotropic thermal displacement tensor ``Uij``.
            j is from (1, 2, 3).

        Returns
        -------
        float
            The elements of anisotropic thermal displacement tensor ``Uij``
        """
        if self.anisotropy:
            return self._U[i, j]
        lat = self.lattice or cartesian_lattice
        return self._U[0, 0] * lat.isotropicunit[i, j]

    def _set_Uij(self, i, j, value):
        """
        Set the elements of anisotropic thermal displacement tensor ``Uij``.

        Parameters
        ----------
        i : int
            The index of anisotropic thermal displacement tensor ``Uij``.
            i is from (1, 2, 3).
        j : int
            The index of anisotropic thermal displacement tensor ``Uij``.
            j is from (1, 2, 3).
        value: float
            The element of anisotropic thermal displacement tensor to be set.
        """
        self._U[i, j] = value
        self._U[j, i] = value
        if not self._anisotropy and i == j != 0:
            self._U[0, 0] = value
        return

    U11 = property(lambda self: self._get_Uij(0, 0),
            lambda self, value: self._set_Uij(0, 0, value), doc =
            "U11 element of anisotropic displacement tensor")
    U22 = property(lambda self: self._get_Uij(1, 1),
            lambda self, value: self._set_Uij(1, 1, value), doc =
            "U22 element of anisotropic displacement tensor")
    U33 = property(lambda self: self._get_Uij(2, 2),
            lambda self, value: self._set_Uij(2, 2, value), doc =
            "U33 element of anisotropic displacement tensor")
    U12 = property(lambda self: self._get_Uij(0, 1),
            lambda self, value: self._set_Uij(0, 1, value), doc =
            "U12 element of anisotropic displacement tensor")
    U13 = property(lambda self: self._get_Uij(0, 2),
            lambda self, value: self._set_Uij(0, 2, value), doc =
            "U13 element of anisotropic displacement tensor")
    U23 = property(lambda self: self._get_Uij(1, 2),
            lambda self, value: self._set_Uij(1, 2, value), doc =
            "U23 element of anisotropic displacement tensor")

    # Uisoequiv

    # TODO: convert to new-style property definition

    @property
    def Uisoequiv(self):
        """
        float: The isotropic displacement parameter.

        When specified anisotropy is set to False.
        Only one of the ``U``, ``Uisoequiv`` arguments can be used.
        """
        if not self.anisotropy:
            return self._U[0, 0]
        if self.lattice is None:
            return numpy.trace(self._U) / 3.0
        lat = self.lattice
        rv = 1.0 / 3.0 * (
                self._U[0,0]*lat.ar*lat.ar*lat.a*lat.a +
                self._U[1,1]*lat.br*lat.br*lat.b*lat.b +
                self._U[2,2]*lat.cr*lat.cr*lat.c*lat.c +
                2*self._U[0,1]*lat.ar*lat.br*lat.a*lat.b*lat.cg +
                2*self._U[0,2]*lat.ar*lat.cr*lat.a*lat.c*lat.cb +
                2*self._U[1,2]*lat.br*lat.cr*lat.b*lat.c*lat.ca)
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

    B11 = property(lambda self: _UtoB*self._get_Uij(0, 0),
            lambda self, value: self._set_Uij(0, 0, _BtoU*value), doc =
            "B11 element of Debye-Waler displacement tensor")
    B22 = property(lambda self: _UtoB*self._get_Uij(1, 1),
            lambda self, value: self._set_Uij(1, 1, _BtoU*value), doc =
            "B22 element of Debye-Waler displacement tensor")
    B33 = property(lambda self: _UtoB*self._get_Uij(2, 2),
            lambda self, value: self._set_Uij(2, 2, _BtoU*value), doc =
            "B33 element of Debye-Waler displacement tensor")
    B12 = property(lambda self: _UtoB*self._get_Uij(0, 1),
            lambda self, value: self._set_Uij(0, 1, _BtoU*value), doc =
            "B12 element of Debye-Waler displacement tensor")
    B13 = property(lambda self: _UtoB*self._get_Uij(0, 2),
            lambda self, value: self._set_Uij(0, 2, _BtoU*value), doc =
            "B13 element of Debye-Waler displacement tensor")
    B23 = property(lambda self: _UtoB*self._get_Uij(1, 2),
            lambda self, value: self._set_Uij(1, 2, _BtoU*value), doc =
            "B23 element of Debye-Waler displacement tensor")

    # Bisoequiv

    # TODO: convert to new-style property definition

    @property
    def Bisoequiv(self):
        """
        float: The Debye-Waler isotropic temperature factor.

        When specified anisotropy is set to False.
        Only one of the ``B``, ``Bisoequiv`` arguments can be used.
        """
        return _UtoB * self.Uisoequiv

    @Bisoequiv.setter
    def Bisoequiv(self, value):
        self.Uisoequiv = _BtoU*value

# End of class Atom

# Local Helpers --------------------------------------------------------------

class _AtomCartesianCoordinates(numpy.ndarray):
    """
    Helper array for accessing Cartesian coordinates.

    Updates of this array are propagated to the xyz fractional coordinates
    of the owner atom according to its lattice attribute.

    Attributes
    ----------
    atom : Atom
        Atom instance linked to these coordinates.

    Parameters
    ----------
    atom : Atom
        Element symbol string or Atom instance.

    """


    def __new__(self, atom):
        """
        Print the atom instance as standard numpy array.

        Returns
        -------
        ndarray
            Printed version of the atom isntance as standard numpy array.
        """
        return numpy.empty(3, dtype=float).view(self)


    def __init__(self, atom):
        """
        Obtain the fractional coordinates within the associated `lattice`, ``xyz``.
        """
        self._atom = atom
        self.asarray[:] = atom.lattice.cartesian(atom.xyz)
        return


    @property
    def asarray(self):
        """
        Print the array as standard numpy array.
        """
        return self.view(numpy.ndarray)


    def __setitem__(self, idx, value):
        """
        Set the idx-th coordinate and update linked self.xyz.

        Parameters
        ----------
        idx : int
            The index in xyz array.
        value: float
            The new value of x, y or z to be set.
        """
        self.asarray[idx] = value
        self._atom.xyz[:] = self._atom.lattice.fractional(self)
        return


    def __array_wrap__(self, out_arr, context=None):
        """
        Let this type of operations yields standard numpy array.

        Parameters
        ----------
        out_arr : ndarray
            Standard numpy array.
        context : bool, optional
            no definition.

        Returns
        -------
        ndarray
            Printed version of the array as standard numpy array.
        """
        return out_arr.view(numpy.ndarray)

# End of _AtomCartesianCoordinates

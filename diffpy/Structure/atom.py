##############################################################################
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
##############################################################################

"""class Atom for storing properties of a single atom"""

__id__ = "$Id$"

import numpy
from diffpy.Structure.lattice import cartesian as cartesian_lattice
from diffpy.Structure import IsotropyError

# conversion constants
_BtoU = 1.0/(8 * numpy.pi**2)
_UtoB = 1.0/_BtoU

class CartesianCoordinatesArray(numpy.ndarray):
    """Helper array for accessing Cartesian coordinates.
    Converts and updates related array of corresponding fractional
    coordinates.

    Data members:
        lattice -- instance of Lattice defining fractional coordinates
        xyz     -- instance of numpy.array storing fractional coordinates
    """

    def __new__(self, lattice, xyz):
        return numpy.zeros(3, dtype=float).view(self)

    def __init__(self, lattice, xyz):
        self.lattice = lattice
        self.xyz = xyz
        self[:] = self.lattice.cartesian(self.xyz)
        pass

    def __setitem__(self, idx, value):
        """Set idx-th coordinate and update linked self.xyz

        idx     -- index in xyz array
        value   -- new value of x, y or z
        """
        numpy.ndarray.__setitem__(self, idx, value)
        self.xyz[:] = self.lattice.fractional(self)
        return

# End of CartesianCoordinatesArray


class Atom(object):
    """Atom --> class for storing atom information

    Data members:
        element     -- type of the atom
        xyz         -- fractional coordinates
        name        -- atom label
        occupancy   -- fractional occupancy
        xyz_cartn   -- absolute Cartesian coordinates, property synced with xyz
        anisotropy  -- flag for anisotropic thermal displacements, property
        U           -- anisotropic thermal displacement tensor, property
        Uij         -- elements of U tensor, where i, j are from (1, 2, 3),
                       property
        Uisoequiv   -- isotropic thermal displacement or equivalent value,
                       property
        Bisoequiv   -- Debye-Waler isotropic temperature factor or equivalent
                       value, property
        lattice     -- coordinate system for fractional coordinates,
                       an instance of Lattice or None for Cartesian system

    Private data:
        _U          -- storage of U property, 3x3 numpy matrix
        _Uisoequiv  -- storage of Uisoequiv property, float
        _anisotropy -- storage of anisotropy property, bool or
                       None when not determined
        _Uijsynced  -- flag for consistency of _U with _Uisoequiv,
                       it is meaningful only for isotropic atoms

    Class data:
        tol_anisotropy -- tolerated U matrix deviation for isotropic atom
    """

    tol_anisotropy = 1.0e-6

    def __init__(self, atype=None, xyz=None, name=None, occupancy=None,
            anisotropy=None, U=None, Uisoequiv=None, lattice=None):
        """Create atom of a specified type at given lattice coordinates.
        Atom(a) creates a copy of Atom instance a.

        atype       -- element symbol string or Atom instance
        xyz         -- fractional coordinates
        name        -- atom label
        occupancy   -- fractional occupancy
        anisotropy  -- flag for anisotropic thermal displacements
        U           -- anisotropic thermal displacement tensor, property
        Uisoequiv   -- isotropic thermal displacement or equivalent value,
                       property
        lattice     -- coordinate system for fractional coordinates
        """
        # declare data members
        self.element = None
        self.xyz = numpy.zeros(3, dtype=float)
        self.name = ''
        self.occupancy = 1.0
        self._anisotropy = None
        self._U = numpy.zeros((3,3), dtype=float)
        self._Uisoequiv = 0.0
        self._Usynced = True
        self.lattice = None
        # assign them as needed
        if isinstance(atype, Atom):
            atype_dup = atype.__copy__()
            self.__dict__.update(atype_dup.__dict__)
        else:
            self.element = atype
        # take care of remaining arguments
        if xyz is not None:         self.xyz[:] = xyz
        if name is not None:        self.name = name
        if occupancy is not None:   self.occupancy = float(occupancy)
        if anisotropy is not None:  self._anisotropy = bool(anisotropy)
        if U is not None:           self._U = U
        if Uisoequiv is not None:   self._Uisoequiv = Uisoequiv
        if lattice is not None:     self.lattice = lattice
        return

    def msdLat(self, vl):
        """mean square displacement of an atom along lattice vector

        vl -- vector in lattice coordinates

        return mean square displacement
        """
        if not self.anisotropy:     return self._Uisoequiv
        # here we need to calculate msd
        lat = self.lattice or cartesian_lattice
        vln = numpy.array(vl, dtype=float)/lat.norm(vl)
        G = lat.metrics
        rhs = numpy.array([ G[0]*lat.ar,
                            G[1]*lat.br,
                            G[2]*lat.cr ], dtype=float)
        rhs = numpy.dot(rhs, vln)
        msd = numpy.dot(rhs, numpy.dot(self._U, rhs))
        return msd

    def msdCart(self, vc):
        """mean square displacement of an atom along cartesian vector

        vc -- vector in absolute cartesian coordinates

        return mean square displacement
        """
        if not self.anisotropy:     return self._Uisoequiv
        # here we need to calculate msd
        lat = self.lattice or cartesian_lattice
        vcn = numpy.array(vc, dtype=float)
        vcn /= numpy.sqrt(numpy.sum(vcn**2))
        F1 = lat.normbase
        Uc = numpy.dot(numpy.transpose(F1), numpy.dot(self._U, F1))
        msd = numpy.dot(vcn, numpy.dot(Uc, vcn))
        return msd

    def __repr__(self):
        """simple string representation"""
        xyz = self.xyz
        s = "%-4s %8.6f %8.6f %8.6f %6.4f" % \
                (self.element, xyz[0], xyz[1], xyz[2], self.occupancy)
        return s

    def __copy__(self):
        """Return a copy of this instance.
        """
        adup = Atom(self.element)
        adup.__dict__.update(self.__dict__)
        # create copies for what should be copied
        adup.xyz = numpy.array(self.xyz)
        adup._U = numpy.array(self._U)
        return adup

    ####################################################################
    # property handlers
    ####################################################################

    # xyz_cartn

    def _get_xyz_cartn(self):
        if not self.lattice:
            rv = self.xyz
        else:
            rv = CartesianCoordinatesArray(self.lattice, self.xyz)
        return rv

    def _set_xyz_cartn(self, value):
        if not self.lattice:
            self.xyz[:] = value
        else:
            self.xyz = self.lattice.fractional(value)
        return

    xyz_cartn = property(_get_xyz_cartn, _set_xyz_cartn, doc =
        """absolute Cartesian coordinates of an atom
        """ )

    # anisotropy

    def _get_anisotropy(self):
        # determine when unknown
        if self._anisotropy is None:
            Uisoequiv = self._get_Uisoequiv()
            # calculate isotropic tensor Uisoij
            lat = self.lattice or cartesian_lattice
            Tu = lat.recnormbase
            Uisoij = numpy.dot(numpy.transpose(Tu), Uisoequiv*Tu)
            # compare with new value
            maxUdiff = numpy.max(numpy.fabs(self._U - Uisoij))
            self._anisotropy = bool(maxUdiff > Atom.tol_anisotropy)
            self._Uijsynced = False
        return self._anisotropy

    def _set_anisotropy(self, value):
        if bool(value) is self._anisotropy: return
        # convert from isotropic to anisotropic
        if value:
            self._U = self._get_U()
        # otherwise convert from anisotropic to isotropic
        else:
            self._Uisoequiv = self._get_Uisoequiv()
            self._Uijsynced = False
        self._anisotropy = bool(value)
        return

    anisotropy = property(_get_anisotropy, _set_anisotropy, doc =
        """flag for anisotropic thermal displacements.
        """ )

    # U

    def _get_U(self):
        # for isotropic non-synced case we need to
        # calculate _U from _Uisoequiv
        if self._anisotropy is False and not self._Uijsynced:
            lat = self.lattice or cartesian_lattice
            Tu = lat.recnormbase
            self._U = numpy.dot(numpy.transpose(Tu), self._Uisoequiv*Tu)
            self._Uijsynced = True
        # handle can be changed by the caller
        self._anisotropy = None
        return self._U

    def _set_U(self, value):
        self._anisotropy = None
        self._U = numpy.array(value, dtype=float)
        return

    U = property(_get_U, _set_U, doc =
        "anisotropic thermal displacement tensor.")

    # Uij elements

    def _get_Uij(self, i, j):
        Uij = self._get_U()
        return Uij[i,j]

    def _set_Uij(self, i, j, value):
        self._anisotropy = None
        self._U[i,j] = value
        self._U[j,i] = value

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

    def _get_Uisoequiv(self):
        if self._anisotropy is None or self._anisotropy is True:
            lat = self.lattice or cartesian_lattice
            Uequiv = (
                    self._U[0,0]*lat.ar*lat.ar*lat.a*lat.a +
                    self._U[1,1]*lat.br*lat.br*lat.b*lat.b +
                    self._U[2,2]*lat.cr*lat.cr*lat.c*lat.c +
                    2*self._U[0,1]*lat.ar*lat.br*lat.a*lat.b*lat.cg +
                    2*self._U[0,2]*lat.ar*lat.cr*lat.a*lat.c*lat.cb +
                    2*self._U[1,2]*lat.br*lat.cr*lat.b*lat.c*lat.ca ) / 3.0
            self._Uisoequiv = Uequiv
        else:
            self._Uisoequiv = self._U[0,0]
        return self._Uisoequiv

    def _set_Uisoequiv(self, value):
        double_eps = (1.0 + numpy.sqrt(2.0**-52)) - 1.0
        self._Uisoequiv = float(value)
        self._Uijsynced = False
        if self._get_anisotropy():
            Uequiv = self._get_Uisoequiv()
            # scale if Uequiv is not zero
            if numpy.fabs(Uequiv) > double_eps:
                self._U *= value/Uequiv
        # otherwise just convert from Uiso value
        else:
            lat = self.lattice or cartesian_lattice
            Tu = lat.recnormbase
            self._U = numpy.dot(numpy.transpose(Tu), value*Tu)
            self._Uijsynced = True
        return

    Uisoequiv = property(_get_Uisoequiv, _set_Uisoequiv, doc =
            "isotropic thermal displacement or equivalent value")

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

    def _get_Bisoequiv(self):
        return _UtoB * self._get_Uisoequiv()

    def _set_Bisoequiv(self, value):
        self._set_Uisoequiv(_BtoU*value)

    Bisoequiv = property(_get_Bisoequiv, _set_Bisoequiv, doc =
            "Debye-Waler isotropic thermal displacement or equivalent value")

# End of class Atom

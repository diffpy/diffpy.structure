#!/usr/bin/env python
# Copyright 2002 by PyMMLib Development Group, http://pymmlib.sourceforge.net/
# This code is part of the PyMMLib distribution and governed by
# its license. Please see the LICENSE_pymmlib file that should have been
# included as part of this package.
"""Symmetry operations as functions on vectors or arrays.
"""

import numpy

# 64 unique rotation matricies
Rot_Z_mY_X = numpy.array([[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]], float)
Rot_Y_mX_mZ = numpy.array([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], float)
Rot_XmY_X_mZ = numpy.array([[1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], float)
Rot_mX_Y_mZ = numpy.array([[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], float)
Rot_X_mZ_Y = numpy.array([[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], float)
Rot_Y_mXY_Z = numpy.array([[0.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], float)
Rot_Y_mX_Z = numpy.array([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], float)
Rot_XmY_X_Z = numpy.array([[1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], float)
Rot_mX_mXY_mZ = numpy.array([[-1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, -1.0]], float)
Rot_Y_Z_X = numpy.array([[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], float)
Rot_mY_mZ_X = numpy.array([[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], float)
Rot_X_Z_mY = numpy.array([[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]], float)
Rot_XmY_mY_Z = numpy.array([[1.0, -1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], float)
Rot_Y_X_mZ = numpy.array([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], float)
Rot_Y_mZ_X = numpy.array([[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]], float)
Rot_mXY_Y_Z = numpy.array([[-1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], float)
Rot_mX_mY_mZ = numpy.array([[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], float)
Rot_X_Y_mZ = numpy.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], float)
Rot_mXY_mX_Z = numpy.array([[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], float)
Rot_mZ_mY_mX = numpy.array([[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [-1.0, 0.0, 0.0]], float)
Rot_X_mZ_mY = numpy.array([[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]], float)
Rot_X_Y_Z = numpy.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], float)
Rot_mY_mX_mZ = numpy.array([[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], float)
Rot_mY_X_Z = numpy.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], float)
Rot_Z_X_Y = numpy.array([[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], float)
Rot_X_XmY_Z = numpy.array([[1.0, 0.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], float)
Rot_mY_X_mZ = numpy.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], float)
Rot_mY_Z_mX = numpy.array([[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]], float)
Rot_mY_Z_X = numpy.array([[0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], float)
Rot_mX_mZ_mY = numpy.array([[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0]], float)
Rot_mX_Z_Y = numpy.array([[-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], float)
Rot_mZ_mX_mY = numpy.array([[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], float)
Rot_X_XmY_mZ = numpy.array([[1.0, 0.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, -1.0]], float)
Rot_mY_XmY_mZ = numpy.array([[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, -1.0]], float)
Rot_Z_X_mY = numpy.array([[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], float)
Rot_mZ_mY_X = numpy.array([[0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [1.0, 0.0, 0.0]], float)
Rot_X_Z_Y = numpy.array([[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]], float)
Rot_Z_mX_mY = numpy.array([[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], float)
Rot_mX_Z_mY = numpy.array([[-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]], float)
Rot_X_mY_Z = numpy.array([[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], float)
Rot_mY_mX_Z = numpy.array([[0.0, -1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], float)
Rot_Z_mY_mX = numpy.array([[0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [-1.0, 0.0, 0.0]], float)
Rot_mX_mY_Z = numpy.array([[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]], float)
Rot_Z_Y_X = numpy.array([[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], float)
Rot_mZ_Y_mX = numpy.array([[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], float)
Rot_Y_Z_mX = numpy.array([[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]], float)
Rot_mY_XmY_Z = numpy.array([[0.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, 1.0]], float)
Rot_mXY_Y_mZ = numpy.array([[-1.0, 1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]], float)
Rot_mZ_mX_Y = numpy.array([[0.0, 0.0, -1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], float)
Rot_mX_mZ_Y = numpy.array([[-1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]], float)
Rot_mX_Y_Z = numpy.array([[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], float)
Rot_X_mY_mZ = numpy.array([[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], float)
Rot_mZ_X_Y = numpy.array([[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], float)
Rot_Y_mZ_mX = numpy.array([[0.0, 1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], float)
Rot_mY_mZ_mX = numpy.array([[0.0, -1.0, 0.0], [0.0, 0.0, -1.0], [-1.0, 0.0, 0.0]], float)
Rot_mZ_Y_X = numpy.array([[0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]], float)
Rot_Z_Y_mX = numpy.array([[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], float)
Rot_mXY_mX_mZ = numpy.array([[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0]], float)
Rot_XmY_mY_mZ = numpy.array([[1.0, -1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], float)
Rot_Z_mX_Y = numpy.array([[0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], float)
Rot_mX_mXY_Z = numpy.array([[-1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 1.0]], float)
Rot_Y_mXY_mZ = numpy.array([[0.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, -1.0]], float)
Rot_mZ_X_mY = numpy.array([[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0]], float)
Rot_Y_X_Z = numpy.array([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], float)

# 32 unique translation vectors
Tr_0_0_34 = numpy.array([0.0, 0.0, 3.0 / 4.0], float)
Tr_12_0_34 = numpy.array([1.0 / 2.0, 0.0, 3.0 / 4.0], float)
Tr_0_0_56 = numpy.array([0.0, 0.0, 5.0 / 6.0], float)
Tr_12_0_12 = numpy.array([1.0 / 2.0, 0.0, 1.0 / 2.0], float)
Tr_0_12_12 = numpy.array([0.0, 1.0 / 2.0, 1.0 / 2.0], float)
Tr_12_0_14 = numpy.array([1.0 / 2.0, 0.0, 1.0 / 4.0], float)
Tr_0_12_14 = numpy.array([0.0, 1.0 / 2.0, 1.0 / 4.0], float)
Tr_14_14_14 = numpy.array([1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0], float)
Tr_0_12_34 = numpy.array([0.0, 1.0 / 2.0, 3.0 / 4.0], float)
Tr_34_14_14 = numpy.array([3.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0], float)
Tr_0_0_0 = numpy.array([0.0, 0.0, 0.0], float)
Tr_23_13_56 = numpy.array([2.0 / 3.0, 1.0 / 3.0, 5.0 / 6.0], float)
Tr_14_14_34 = numpy.array([1.0 / 4.0, 1.0 / 4.0, 3.0 / 4.0], float)
Tr_12_12_0 = numpy.array([1.0 / 2.0, 1.0 / 2.0, 0.0], float)
Tr_23_13_13 = numpy.array([2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0], float)
Tr_13_23_23 = numpy.array([1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0], float)
Tr_12_12_12 = numpy.array([1.0 / 2.0, 1.0 / 2.0, 1.0 / 2.0], float)
Tr_12_12_14 = numpy.array([1.0 / 2.0, 1.0 / 2.0, 1.0 / 4.0], float)
Tr_14_34_14 = numpy.array([1.0 / 4.0, 3.0 / 4.0, 1.0 / 4.0], float)
Tr_12_12_34 = numpy.array([1.0 / 2.0, 1.0 / 2.0, 3.0 / 4.0], float)
Tr_0_0_23 = numpy.array([0.0, 0.0, 2.0 / 3.0], float)
Tr_0_12_0 = numpy.array([0.0, 1.0 / 2.0, 0.0], float)
Tr_14_34_34 = numpy.array([1.0 / 4.0, 3.0 / 4.0, 3.0 / 4.0], float)
Tr_34_34_14 = numpy.array([3.0 / 4.0, 3.0 / 4.0, 1.0 / 4.0], float)
Tr_12_0_0 = numpy.array([1.0 / 2.0, 0.0, 0.0], float)
Tr_34_34_34 = numpy.array([3.0 / 4.0, 3.0 / 4.0, 3.0 / 4.0], float)
Tr_0_0_13 = numpy.array([0.0, 0.0, 1.0 / 3.0], float)
Tr_0_0_12 = numpy.array([0.0, 0.0, 1.0 / 2.0], float)
Tr_13_23_16 = numpy.array([1.0 / 3.0, 2.0 / 3.0, 1.0 / 6.0], float)
Tr_0_0_14 = numpy.array([0.0, 0.0, 1.0 / 4.0], float)
Tr_0_0_16 = numpy.array([0.0, 0.0, 1.0 / 6.0], float)
Tr_34_14_34 = numpy.array([3.0 / 4.0, 1.0 / 4.0, 3.0 / 4.0], float)


class SymOp(object):
    """The transformation of coordinates to a symmetry-related position.

    The SymOp operation involves rotation and translation in cell coordinates.

    Parameters
    ----------
    R : numpy.ndarray
        The 3x3 matrix of rotation for this symmetry operation.
    t : numpy.ndarray
        The vector of translation in this symmetry operation.

    Attributes
    ----------
    R : numpy.ndarray
        The 3x3 matrix of rotation pertaining to unit cell coordinates.
        This may be identity, simple rotation, improper rotation, mirror
        or inversion. The determinant of *R* is either +1 or -1.
    t : numpy.ndarray
        The translation of cell coordinates applied after rotation *R*.
    """

    def __init__(self, R, t):
        self.R = R
        self.t = t
        return

    def __str__(self):
        """Printable representation of this SymOp object."""
        x = "[%6.3f %6.3f %6.3f %6.3f]\n" % (self.R[0, 0], self.R[0, 1], self.R[0, 2], self.t[0])
        x += "[%6.3f %6.3f %6.3f %6.3f]\n" % (self.R[1, 0], self.R[1, 1], self.R[1, 2], self.t[1])
        x += "[%6.3f %6.3f %6.3f %6.3f]\n" % (self.R[2, 0], self.R[2, 1], self.R[2, 2], self.t[2])
        return x

    def __call__(self, vec):
        """Return symmetry-related position for the specified coordinates.

        Parameters
        ----------
        vec : numpy.ndarray
            The initial position in fractional cell coordinates.

        Returns
        -------
        numpy.ndarray
            The transformed position after this symmetry operation.
        """
        return numpy.dot(self.R, vec) + self.t

    def __eq__(self, symop):
        """Implement the ``self == symop`` test of equality.

        Return ``True`` when *self* and *symop* difference is within
        tiny round-off errors.
        """
        return numpy.allclose(self.R, symop.R) and numpy.allclose(self.t, symop.t)

    def is_identity(self):
        """Check if this SymOp is an identity operation.

        Returns
        -------
        bool
            ``True`` if this is an identity operation within a small round-off.
            Return ``False`` otherwise.
        """
        rv = numpy.allclose(self.R, numpy.identity(3, float)) and numpy.allclose(self.t, numpy.zeros(3, float))
        return rv


# End of class SymOp


class SpaceGroup(object):
    """Definition and basic operations for a specific space group.

    Provide standard names and all symmetry operations contained in
    one space group.

    Parameters
    ----------
    number : int
        The space group number.
    num_sym_equiv : int
        The number of symmetry equivalent sites for a general position.
    num_primitive_sym_equiv : int
        The number of symmetry equivalent sites in a primitive unit cell.
    short_name : str
        The short Hermann-Mauguin symbol of the space group.
    point_group_name : str
        The point group of this space group.
    crystal_system : str
        The crystal system of this space group.
    pdb_name : str
        The full Hermann-Mauguin symbol of the space group.
    symop_list : list of SymOp
        The symmetry operations contained in this space group.

    Attributes
    ----------
    number : int
        A unique space group number. This may be incremented by
        several thousands to facilitate unique values for multiple
        settings of the same space group. Use ``number % 1000``
        to get the standard space group number from International
        Tables.
    num_sym_equiv : int
        The number of symmetry equivalent sites for a general position.
    num_primitive_sym_equiv : int
        The number of symmetry equivalent sites in a primitive unit cell.
    short_name : str
        The short Hermann-Mauguin symbol of the space group.
    point_group_name : str
        The point group to which this space group belongs to.
    crystal_system : str
        The crystal system of this space group. The possible values are
        ``"TRICLINIC", "MONOCLINIC", "ORTHORHOMBIC", "TETRAGONAL",
        "TRIGONAL" "HEXAGONAL", "CUBIC"``.
    pdb_name : str
        The full Hermann-Mauguin symbol of the space group.
    symop_list : list of SymOp
        A list of `SymOp` objects for all symmetry operations
        in this space group.
    """

    def __init__(
        self,
        number=None,
        num_sym_equiv=None,
        num_primitive_sym_equiv=None,
        short_name=None,
        point_group_name=None,
        crystal_system=None,
        pdb_name=None,
        symop_list=None,
    ):

        self.number = number
        self.num_sym_equiv = num_sym_equiv
        self.num_primitive_sym_equiv = num_primitive_sym_equiv
        self.short_name = short_name
        self.point_group_name = point_group_name
        self.crystal_system = crystal_system
        self.pdb_name = pdb_name
        self.symop_list = symop_list

    def __repr__(self):
        """Return a string representation of the space group."""
        return ("SpaceGroup #%i (%s, %s). Symmetry matrices: %i, " "point sym. matr.: %i") % (
            self.number,
            self.short_name,
            self.crystal_system[0] + self.crystal_system[1:].lower(),
            self.num_sym_equiv,
            self.num_primitive_sym_equiv,
        )

    def iter_symops(self):
        """Iterate over all symmetry operations in the space group.

        Yields
        ------
        SymOp
            Generate all symmetry operations for this space group.
        """
        return iter(self.symop_list)

    def check_group_name(self, name):
        """Check if given name matches this space group.

        Parameters
        ----------
        name : str or int
            The space group identifier, a string name or number.

        Returns
        -------
        bool
            ``True`` if the specified name matches one of the recognized
            names of this space group or if it equals its `number`.
            Return ``False`` otherwise.
        """

        if name == self.short_name:
            return True
        if name == self.pdb_name:
            return True
        if name == self.point_group_name:
            return True
        if name == self.number:
            return True
        return False

    def iter_equivalent_positions(self, vec):
        """Generate symmetry equivalent positions for the specified position.

        The initial position must be in fractional coordinates and so
        are the symmetry equivalent positions yielded by iteration.
        This generates `num_sym_equiv` positions regardless of initial
        coordinates being a special symmetry position or not.

        Parameters
        ----------
        vec : numpy.ndarray
            The initial position in fractional coordinates.

        Yields
        ------
        numpy.ndarray
            The symmetry equivalent positions in fractional coordinates.
            The positions may be duplicate or outside of the ``0 <= x < 1``
            unit cell bounds.
        """
        for symop in self.symop_list:
            yield symop(vec)
        pass


# End of class SpaceGroup

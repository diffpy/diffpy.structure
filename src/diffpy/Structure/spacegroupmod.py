#!/usr/bin/env python
## Copyright 2002 by PyMMLib Development Group, http://pymmlib.sourceforge.net/
## This code is part of the PyMMLib distribution and governed by
## its license.  Please see the LICENSE_pymmlib file that should have been
## included as part of this package.
"""Symmetry operations as functions on vectors or arrays.
"""

import numpy

## 64 unique rotation matricies
Rot_Z_mY_X    = numpy.array([[ 0.0, 0.0, 1.0], [ 0.0,-1.0, 0.0], [ 1.0, 0.0, 0.0]], float)
Rot_Y_mX_mZ   = numpy.array([[ 0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [ 0.0, 0.0,-1.0]], float)
Rot_XmY_X_mZ  = numpy.array([[ 1.0,-1.0, 0.0], [ 1.0, 0.0, 0.0], [ 0.0, 0.0,-1.0]], float)
Rot_mX_Y_mZ   = numpy.array([[-1.0, 0.0, 0.0], [ 0.0, 1.0, 0.0], [ 0.0, 0.0,-1.0]], float)
Rot_X_mZ_Y    = numpy.array([[ 1.0, 0.0, 0.0], [ 0.0, 0.0,-1.0], [ 0.0, 1.0, 0.0]], float)
Rot_Y_mXY_Z   = numpy.array([[ 0.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [ 0.0, 0.0, 1.0]], float)
Rot_Y_mX_Z    = numpy.array([[ 0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [ 0.0, 0.0, 1.0]], float)
Rot_XmY_X_Z   = numpy.array([[ 1.0,-1.0, 0.0], [ 1.0, 0.0, 0.0], [ 0.0, 0.0, 1.0]], float)
Rot_mX_mXY_mZ = numpy.array([[-1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [ 0.0, 0.0,-1.0]], float)
Rot_Y_Z_X     = numpy.array([[ 0.0, 1.0, 0.0], [ 0.0, 0.0, 1.0], [ 1.0, 0.0, 0.0]], float)
Rot_mY_mZ_X   = numpy.array([[ 0.0,-1.0, 0.0], [ 0.0, 0.0,-1.0], [ 1.0, 0.0, 0.0]], float)
Rot_X_Z_mY    = numpy.array([[ 1.0, 0.0, 0.0], [ 0.0, 0.0, 1.0], [ 0.0,-1.0, 0.0]], float)
Rot_XmY_mY_Z  = numpy.array([[ 1.0,-1.0, 0.0], [ 0.0,-1.0, 0.0], [ 0.0, 0.0, 1.0]], float)
Rot_Y_X_mZ    = numpy.array([[ 0.0, 1.0, 0.0], [ 1.0, 0.0, 0.0], [ 0.0, 0.0,-1.0]], float)
Rot_Y_mZ_X    = numpy.array([[ 0.0, 1.0, 0.0], [ 0.0, 0.0,-1.0], [ 1.0, 0.0, 0.0]], float)
Rot_mXY_Y_Z   = numpy.array([[-1.0, 1.0, 0.0], [ 0.0, 1.0, 0.0], [ 0.0, 0.0, 1.0]], float)
Rot_mX_mY_mZ  = numpy.array([[-1.0, 0.0, 0.0], [ 0.0,-1.0, 0.0], [ 0.0, 0.0,-1.0]], float)
Rot_X_Y_mZ    = numpy.array([[ 1.0, 0.0, 0.0], [ 0.0, 1.0, 0.0], [ 0.0, 0.0,-1.0]], float)
Rot_mXY_mX_Z  = numpy.array([[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [ 0.0, 0.0, 1.0]], float)
Rot_mZ_mY_mX  = numpy.array([[ 0.0, 0.0,-1.0], [ 0.0,-1.0, 0.0], [-1.0, 0.0, 0.0]], float)
Rot_X_mZ_mY   = numpy.array([[ 1.0, 0.0, 0.0], [ 0.0, 0.0,-1.0], [ 0.0,-1.0, 0.0]], float)
Rot_X_Y_Z     = numpy.array([[ 1.0, 0.0, 0.0], [ 0.0, 1.0, 0.0], [ 0.0, 0.0, 1.0]], float)
Rot_mY_mX_mZ  = numpy.array([[ 0.0,-1.0, 0.0], [-1.0, 0.0, 0.0], [ 0.0, 0.0,-1.0]], float)
Rot_mY_X_Z    = numpy.array([[ 0.0,-1.0, 0.0], [ 1.0, 0.0, 0.0], [ 0.0, 0.0, 1.0]], float)
Rot_Z_X_Y     = numpy.array([[ 0.0, 0.0, 1.0], [ 1.0, 0.0, 0.0], [ 0.0, 1.0, 0.0]], float)
Rot_X_XmY_Z   = numpy.array([[ 1.0, 0.0, 0.0], [ 1.0,-1.0, 0.0], [ 0.0, 0.0, 1.0]], float)
Rot_mY_X_mZ   = numpy.array([[ 0.0,-1.0, 0.0], [ 1.0, 0.0, 0.0], [ 0.0, 0.0,-1.0]], float)
Rot_mY_Z_mX   = numpy.array([[ 0.0,-1.0, 0.0], [ 0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]], float)
Rot_mY_Z_X    = numpy.array([[ 0.0,-1.0, 0.0], [ 0.0, 0.0, 1.0], [ 1.0, 0.0, 0.0]], float)
Rot_mX_mZ_mY  = numpy.array([[-1.0, 0.0, 0.0], [ 0.0, 0.0,-1.0], [ 0.0,-1.0, 0.0]], float)
Rot_mX_Z_Y    = numpy.array([[-1.0, 0.0, 0.0], [ 0.0, 0.0, 1.0], [ 0.0, 1.0, 0.0]], float)
Rot_mZ_mX_mY  = numpy.array([[ 0.0, 0.0,-1.0], [-1.0, 0.0, 0.0], [ 0.0,-1.0, 0.0]], float)
Rot_X_XmY_mZ  = numpy.array([[ 1.0, 0.0, 0.0], [ 1.0,-1.0, 0.0], [ 0.0, 0.0,-1.0]], float)
Rot_mY_XmY_mZ = numpy.array([[ 0.0,-1.0, 0.0], [ 1.0,-1.0, 0.0], [ 0.0, 0.0,-1.0]], float)
Rot_Z_X_mY    = numpy.array([[ 0.0, 0.0, 1.0], [ 1.0, 0.0, 0.0], [ 0.0,-1.0, 0.0]], float)
Rot_mZ_mY_X   = numpy.array([[ 0.0, 0.0,-1.0], [ 0.0,-1.0, 0.0], [ 1.0, 0.0, 0.0]], float)
Rot_X_Z_Y     = numpy.array([[ 1.0, 0.0, 0.0], [ 0.0, 0.0, 1.0], [ 0.0, 1.0, 0.0]], float)
Rot_Z_mX_mY   = numpy.array([[ 0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [ 0.0,-1.0, 0.0]], float)
Rot_mX_Z_mY   = numpy.array([[-1.0, 0.0, 0.0], [ 0.0, 0.0, 1.0], [ 0.0,-1.0, 0.0]], float)
Rot_X_mY_Z    = numpy.array([[ 1.0, 0.0, 0.0], [ 0.0,-1.0, 0.0], [ 0.0, 0.0, 1.0]], float)
Rot_mY_mX_Z   = numpy.array([[ 0.0,-1.0, 0.0], [-1.0, 0.0, 0.0], [ 0.0, 0.0, 1.0]], float)
Rot_Z_mY_mX   = numpy.array([[ 0.0, 0.0, 1.0], [ 0.0,-1.0, 0.0], [-1.0, 0.0, 0.0]], float)
Rot_mX_mY_Z   = numpy.array([[-1.0, 0.0, 0.0], [ 0.0,-1.0, 0.0], [ 0.0, 0.0, 1.0]], float)
Rot_Z_Y_X     = numpy.array([[ 0.0, 0.0, 1.0], [ 0.0, 1.0, 0.0], [ 1.0, 0.0, 0.0]], float)
Rot_mZ_Y_mX   = numpy.array([[ 0.0, 0.0,-1.0], [ 0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], float)
Rot_Y_Z_mX    = numpy.array([[ 0.0, 1.0, 0.0], [ 0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]], float)
Rot_mY_XmY_Z  = numpy.array([[ 0.0,-1.0, 0.0], [ 1.0,-1.0, 0.0], [ 0.0, 0.0, 1.0]], float)
Rot_mXY_Y_mZ  = numpy.array([[-1.0, 1.0, 0.0], [ 0.0, 1.0, 0.0], [ 0.0, 0.0,-1.0]], float)
Rot_mZ_mX_Y   = numpy.array([[ 0.0, 0.0,-1.0], [-1.0, 0.0, 0.0], [ 0.0, 1.0, 0.0]], float)
Rot_mX_mZ_Y   = numpy.array([[-1.0, 0.0, 0.0], [ 0.0, 0.0,-1.0], [ 0.0, 1.0, 0.0]], float)
Rot_mX_Y_Z    = numpy.array([[-1.0, 0.0, 0.0], [ 0.0, 1.0, 0.0], [ 0.0, 0.0, 1.0]], float)
Rot_X_mY_mZ   = numpy.array([[ 1.0, 0.0, 0.0], [ 0.0,-1.0, 0.0], [ 0.0, 0.0,-1.0]], float)
Rot_mZ_X_Y    = numpy.array([[ 0.0, 0.0,-1.0], [ 1.0, 0.0, 0.0], [ 0.0, 1.0, 0.0]], float)
Rot_Y_mZ_mX   = numpy.array([[ 0.0, 1.0, 0.0], [ 0.0, 0.0,-1.0], [-1.0, 0.0, 0.0]], float)
Rot_mY_mZ_mX  = numpy.array([[ 0.0,-1.0, 0.0], [ 0.0, 0.0,-1.0], [-1.0, 0.0, 0.0]], float)
Rot_mZ_Y_X    = numpy.array([[ 0.0, 0.0,-1.0], [ 0.0, 1.0, 0.0], [ 1.0, 0.0, 0.0]], float)
Rot_Z_Y_mX    = numpy.array([[ 0.0, 0.0, 1.0], [ 0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]], float)
Rot_mXY_mX_mZ = numpy.array([[-1.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [ 0.0, 0.0,-1.0]], float)
Rot_XmY_mY_mZ = numpy.array([[ 1.0,-1.0, 0.0], [ 0.0,-1.0, 0.0], [ 0.0, 0.0,-1.0]], float)
Rot_Z_mX_Y    = numpy.array([[ 0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [ 0.0, 1.0, 0.0]], float)
Rot_mX_mXY_Z  = numpy.array([[-1.0, 0.0, 0.0], [-1.0, 1.0, 0.0], [ 0.0, 0.0, 1.0]], float)
Rot_Y_mXY_mZ  = numpy.array([[ 0.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [ 0.0, 0.0,-1.0]], float)
Rot_mZ_X_mY   = numpy.array([[ 0.0, 0.0,-1.0], [ 1.0, 0.0, 0.0], [ 0.0,-1.0, 0.0]], float)
Rot_Y_X_Z     = numpy.array([[ 0.0, 1.0, 0.0], [ 1.0, 0.0, 0.0], [ 0.0, 0.0, 1.0]], float)

## 32 unique translation vectors
Tr_0_0_34     = numpy.array([     0.0,     0.0, 3.0/4.0 ], float)
Tr_12_0_34    = numpy.array([ 1.0/2.0,     0.0, 3.0/4.0 ], float)
Tr_0_0_56     = numpy.array([     0.0,     0.0, 5.0/6.0 ], float)
Tr_12_0_12    = numpy.array([ 1.0/2.0,     0.0, 1.0/2.0 ], float)
Tr_0_12_12    = numpy.array([     0.0, 1.0/2.0, 1.0/2.0 ], float)
Tr_12_0_14    = numpy.array([ 1.0/2.0,     0.0, 1.0/4.0 ], float)
Tr_0_12_14    = numpy.array([     0.0, 1.0/2.0, 1.0/4.0 ], float)
Tr_14_14_14   = numpy.array([ 1.0/4.0, 1.0/4.0, 1.0/4.0 ], float)
Tr_0_12_34    = numpy.array([     0.0, 1.0/2.0, 3.0/4.0 ], float)
Tr_34_14_14   = numpy.array([ 3.0/4.0, 1.0/4.0, 1.0/4.0 ], float)
Tr_0_0_0      = numpy.array([     0.0,     0.0,     0.0 ], float)
Tr_23_13_56   = numpy.array([ 2.0/3.0, 1.0/3.0, 5.0/6.0 ], float)
Tr_14_14_34   = numpy.array([ 1.0/4.0, 1.0/4.0, 3.0/4.0 ], float)
Tr_12_12_0    = numpy.array([ 1.0/2.0, 1.0/2.0,     0.0 ], float)
Tr_23_13_13   = numpy.array([ 2.0/3.0, 1.0/3.0, 1.0/3.0 ], float)
Tr_13_23_23   = numpy.array([ 1.0/3.0, 2.0/3.0, 2.0/3.0 ], float)
Tr_12_12_12   = numpy.array([ 1.0/2.0, 1.0/2.0, 1.0/2.0 ], float)
Tr_12_12_14   = numpy.array([ 1.0/2.0, 1.0/2.0, 1.0/4.0 ], float)
Tr_14_34_14   = numpy.array([ 1.0/4.0, 3.0/4.0, 1.0/4.0 ], float)
Tr_12_12_34   = numpy.array([ 1.0/2.0, 1.0/2.0, 3.0/4.0 ], float)
Tr_0_0_23     = numpy.array([     0.0,     0.0, 2.0/3.0 ], float)
Tr_0_12_0     = numpy.array([     0.0, 1.0/2.0,     0.0 ], float)
Tr_14_34_34   = numpy.array([ 1.0/4.0, 3.0/4.0, 3.0/4.0 ], float)
Tr_34_34_14   = numpy.array([ 3.0/4.0, 3.0/4.0, 1.0/4.0 ], float)
Tr_12_0_0     = numpy.array([ 1.0/2.0,     0.0,     0.0 ], float)
Tr_34_34_34   = numpy.array([ 3.0/4.0, 3.0/4.0, 3.0/4.0 ], float)
Tr_0_0_13     = numpy.array([     0.0,     0.0, 1.0/3.0 ], float)
Tr_0_0_12     = numpy.array([     0.0,     0.0, 1.0/2.0 ], float)
Tr_13_23_16   = numpy.array([ 1.0/3.0, 2.0/3.0, 1.0/6.0 ], float)
Tr_0_0_14     = numpy.array([     0.0,     0.0, 1.0/4.0 ], float)
Tr_0_0_16     = numpy.array([     0.0,     0.0, 1.0/6.0 ], float)
Tr_34_14_34   = numpy.array([ 3.0/4.0, 1.0/4.0, 3.0/4.0 ], float)


class SymOp(object):
    """A subclass of the tuple class for performing one symmetry operation.
    """
    def __init__(self, R, t):
        self.R = R
        self.t = t

    def __str__(self):
        x  = "[%6.3f %6.3f %6.3f %6.3f]\n" % (
            self.R[0,0], self.R[0,1], self.R[0,2], self.t[0])
        x += "[%6.3f %6.3f %6.3f %6.3f]\n" % (
            self.R[1,0], self.R[1,1], self.R[1,2], self.t[1])
        x += "[%6.3f %6.3f %6.3f %6.3f]\n" % (
            self.R[2,0], self.R[2,1], self.R[2,2], self.t[2])
        return x

    def __call__(self, vec):
        """Return the symmetry operation on argument vector and.
        """
        return numpy.dot(self.R, vec) + self.t

    def __eq__(self, symop):
        return numpy.allclose(self.R, symop.R) and numpy.allclose(self.t, symop.t)

    def is_identity(self):
        """Returns True if this SymOp is an identity symmetry operation
        (no rotation, no translation), otherwise returns False.
        """
        rv = (numpy.allclose(self.R, numpy.identity(3, float)) and
              numpy.allclose(self.t, numpy.zeros(3, float)))
        return rv

# End of class SymOp


class SpaceGroup(object):
    """Contains the various names and symmetry operations for one space
    group.
    """
    def __init__(self,
                 number                  = None,
                 num_sym_equiv           = None,
                 num_primitive_sym_equiv = None,
                 short_name              = None,
                 alt_name                = None,
                 point_group_name        = None,
                 crystal_system          = None,
                 pdb_name                = None,
                 symop_list              = None):

        self.number                  = number
        self.num_sym_equiv           = num_sym_equiv
        self.num_primitive_sym_equiv = num_primitive_sym_equiv
        self.short_name              = short_name
        self.alt_name                = alt_name
        self.point_group_name        = point_group_name
        self.crystal_system          = crystal_system
        self.pdb_name                = pdb_name
        self.symop_list              = symop_list

    def iter_symops(self):
        """Iterates over all SymOps in the SpaceGroup.
        """
        return iter(self.symop_list)

    def check_group_name(self, name):
        """Checks if the given name is a name for this space group,
        returns True or False.  The space group name can be in several forms:
        the short name, the longer PDB-style name, or the space group number.
        """
        if name == self.short_name:       return True
        if name == self.alt_name:         return True
        if name == self.pdb_name:         return True
        if name == self.point_group_name: return True
        if name == self.number:           return True
        return False

    def iter_equivalent_positions(self, vec):
        """Iterate the symmetry equivalent positions of the argument vector.
        The vector must already be in fractional coordinates, and the
        symmetry equivalent vectors are also in fractional coordinates.
        """
        for symop in self.symop_list:
            yield symop(vec)
        pass

# End of class SpaceGroup

# End of file

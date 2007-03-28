########################################################################
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
########################################################################

"""class Lattice stores properites and provides simple operations in lattice
coordinate system.
"""

__id__ = "$Id$"

import copy
import math
import types
import numpy
import numpy.linalg as numalg
from StructureErrors import InvalidLattice

##############################################################################
# helper functions

def cosd(x):
    """return the cosine of x (measured in degrees)"""
    c = {  0.0 : +1.0,   60.0 : +0.5,   90.0 : 0.0,  120.0 : -0.5,
         180.0 : -1.0,  240.0 : -0.5,  270.0 : 0.0,  300.0 : +0.5 }
    return c.get(x % 360.0, math.cos(math.radians(x)))

def sind(x):
    """return the sine of x (measured in degrees)"""
    return cosd(90.0 - x)

# End of helper functions


##############################################################################
class Lattice:
    """Lattice --> general coordinate system

    Data members:

        a, b, c, alpha, beta, gamma -- read-only lattice parameters,
               their values are set by setLatPar() and setLatBase()
        ar, br, cr, alphar, betar, gammar -- read-only parameters of
               reciprocal lattice, set by setLatPar() and setLatBase()
        metrics  -- metrics tensor
        base     -- matrix of base vectors in cartesian coordinates,
                    base = stdbase*baserot
        stdbase  -- matrix of base vectors in standard orientation
        baserot  -- base rotation matrix
        recbase  -- inverse of base matrix, its columns are reciprocal
                    vectors in cartesian coordinates
        normbase -- base with magnitudes of reciprocal vectors
        recnormbase -- inverse of normbase

    Note: All data members are read-only, their values get set by calling
    setLatPar() or setLatBase() methods.
    """

    def __init__(self, a=None, b=None, c=None,
            alpha=None, beta=None, gamma=None,
            baserot=numpy.identity(3, dtype=float), base=None):
        """define new coordinate system, the default is Cartesian
        There are 4 ways how to create Lattice instance:

        Lattice()         -- create cartesian coordinates
        Lattice(a, b, c, alpha, beta, gamma) -- define coordinate system
                             from specified lattice parameters
        Lattice(base=abc) -- create coordinate system using the given base,
                             abc is a 3x3 matrix (or nested list), of row
                             base vectors
        Lattice(lat)      -- create a copy of existing Lattice lat
        """
        # initialize data members, they values will be set by setLatPar()
        self.a = self.b = self.c = None
        self.alpha = self.beta = self.gamma = None
        self.ar = self.br = self.cr = None
        self.alphar = self.betar = self.gammar = None
        self.baserot = None
        self.base = self.recbase = None
        self.normbase = self.recnormbase = None
        # work out argument variants
        # Lattice()
        if [a,b,c,alpha,beta,gamma,base] == 7*[None]:
            self.setLatPar(1.0, 1.0, 1.0, 90.0, 90.0, 90.0, baserot)
        # Lattice(base=abc)
        elif base is not None:
            self.setLatBase(base)
        # Lattice(lat)
        elif isinstance(a, Lattice):
            copythese = dict.fromkeys( ('baserot', 'base',
                'recbase', 'normbase', 'recnormbase') )
            for attribute, value in a.__dict__.iteritems():
                if attribute in copythese:
                    setattr(self, attribute, copy.deepcopy(value))
                    continue
                setattr(self, attribute, value)
        # otherwise do default Lattice(a, b, c, alpha, beta, gamma)
        else:
            self.setLatPar( float(a), float(b), float(c),
                    float(alpha), float(beta), float(gamma), baserot )
        return

    def setLatPar(self, a=None, b=None, c=None,
            alpha=None, beta=None, gamma=None, baserot=None):
        """set lattice parameters and all related tensors

        a, b, c, alpha, beta, gamma -- lattice parameters
        baserot -- unit cell rotation, base = stdbase*baserot

        Note: parameters with value None will remain unchanged.

        return self
        """
        if a is not None: self.a = float(a)
        if b is not None: self.b = float(b)
        if c is not None: self.c = float(c)
        if alpha is not None: self.alpha = float(alpha)
        if beta is not None: self.beta = float(beta)
        if gamma is not None: self.gamma = float(gamma)
        if baserot is not None: self.baserot = numpy.array(baserot)
        (ca, sa) = ( cosd(self.alpha), sind(self.alpha) )
        (cb, sb) = ( cosd(self.beta),  sind(self.beta) )
        (cg, sg) = ( cosd(self.gamma), sind(self.gamma) )
        # Vunit is a volume of unit cell with a=b=c=1
        Vunit = math.sqrt(1.0 + 2.0*ca*cb*cg - ca*ca - cb*cb - cg*cg)
        # reciprocal lattice
        self.ar = sa/(self.a*Vunit)
        self.br = sb/(self.b*Vunit)
        self.cr = sg/(self.c*Vunit)
        car = (cb*cg - ca)/(sb*sg); sar = math.sqrt(1.0 - car**2)
        cbr = (ca*cg - cb)/(sa*sg); sbr = math.sqrt(1.0 - cbr**2)
        cgr = (ca*cb - cg)/(sa*sb); sgr = math.sqrt(1.0 - cgr**2)
        self.alphar = math.degrees(math.acos(car))
        self.betar = math.degrees(math.acos(cbr))
        self.gammar = math.degrees(math.acos(cgr))
        # metric tensor
        self.metrics = numpy.array( [
                [ self.a*self.a,     self.a*self.b*cg,  self.a*self.c*cb ],
                [ self.b*self.a*cg,  self.b*self.b,     self.b*self.c*ca ],
                [ self.c*self.a*cb,  self.c*self.b*ca,  self.c*self.c    ] ],
                dtype=float )
        # standard cartesian coordinates of lattice vectors
        self.stdbase = numpy.array( [
                [ 1.0/self.ar, -cgr/sgr/self.ar, cb*self.a ],
                [ 0.0,         self.b*sa,        self.b*ca ],
                [ 0.0,         0.0,              self.c    ] ],
                dtype=float )
        # cartesian coordinates of lattice vectors
        self.base = numpy.dot(self.stdbase, self.baserot)
        self.recbase = numalg.inv(self.base)
        # bases normalized to unit reciprocal vectors
        self.normbase = numpy.array([ self.base[0,:]*self.ar,
                                    self.base[1,:]*self.br,
                                    self.base[2,:]*self.cr ])
        self.recnormbase = numpy.array(self.recbase)
        self.recnormbase[:,0] /= self.ar
        self.recnormbase[:,1] /= self.br
        self.recnormbase[:,2] /= self.cr
        return self

    def setLatBase(self, base):
        """set matrix of unit cell base vectors and calculate corresponding
        lattice parameters and stdbase, baserot and metrics tensors

        return self
        """
        self.base = numpy.array(base)
        detbase = numalg.det(self.base)
        if abs(detbase) < 1.0e-8:
            raise InvalidLattice, "base vectors are degenerate"
        elif detbase < 0.0:
            raise InvalidLattice, "base is not right-handed"
        self.a = numpy.sqrt(numpy.dot(self.base[0,:], self.base[0,:]))
        self.b = numpy.sqrt(numpy.dot(self.base[1,:], self.base[1,:]))
        self.c = numpy.sqrt(numpy.dot(self.base[2,:], self.base[2,:]))
        ca = numpy.dot(self.base[1,:], self.base[2,:]) / (self.b*self.c)
        cb = numpy.dot(self.base[0,:], self.base[2,:]) / (self.a*self.c)
        cg = numpy.dot(self.base[0,:], self.base[1,:]) / (self.a*self.b)
        sa = numpy.sqrt(1.0 - ca**2)
        sb = numpy.sqrt(1.0 - cb**2)
        sg = numpy.sqrt(1.0 - cg**2)
        self.alpha = math.degrees(math.acos(ca))
        self.beta = math.degrees(math.acos(cb))
        self.gamma = math.degrees(math.acos(cg))
        # Vunit is a volume of unit cell with a=b=c=1
        Vunit = math.sqrt(1.0 + 2.0*ca*cb*cg - ca*ca - cb*cb - cg*cg)
        # reciprocal lattice
        self.ar = sa/(self.a*Vunit)
        self.br = sb/(self.b*Vunit)
        self.cr = sg/(self.c*Vunit)
        car = (cb*cg - ca)/(sb*sg); sar = math.sqrt(1.0 - car**2)
        cbr = (ca*cg - cb)/(sa*sg); sbr = math.sqrt(1.0 - cbr**2)
        cgr = (ca*cb - cg)/(sa*sb); sgr = math.sqrt(1.0 - cgr**2)
        self.alphar = math.degrees(math.acos(car))
        self.betar = math.degrees(math.acos(cbr))
        self.gammar = math.degrees(math.acos(cgr))
        # standard orientation of lattice vectors
        self.stdbase = numpy.array( [
                [ 1.0/self.ar, -cgr/sgr/self.ar, cb*self.a ],
                [ 0.0,         self.b*sa,        self.b*ca ],
                [ 0.0,         0.0,              self.c    ] ],
                dtype=float )
        # calculate unit cell rotation matrix,  base = stdbase*baserot
        self.baserot = numpy.dot( numalg.inv(self.stdbase), self.base )
        self.recbase = numalg.inv(self.base)
        # bases normalized to unit reciprocal vectors
        self.normbase = numpy.array([ self.base[0,:]*self.ar,
                                    self.base[1,:]*self.br,
                                    self.base[2,:]*self.cr ])
        self.recnormbase = numpy.array(self.recbase)
        self.recnormbase[:,0] /= self.ar
        self.recnormbase[:,1] /= self.br
        self.recnormbase[:,2] /= self.cr
        # update metrics tensor
        self.metrics = numpy.array( [
                [ self.a*self.a,     self.a*self.b*cg,  self.a*self.c*cb ],
                [ self.b*self.a*cg,  self.b*self.b,     self.b*self.c*ca ],
                [ self.c*self.a*cb,  self.c*self.b*ca,  self.c*self.c    ] ],
                dtype=float )
        return self

    def abcABG(self):
        """Return tuple of 6 lattice parameters.
        """
        rv = (self.a, self.b, self.c, self.alpha, self.beta, self.gamma)
        return rv

    def reciprocal(self):
        """return the reciprocal lattice to self"""
        from copy import deepcopy
        rec = deepcopy(self)
        rec.setLatBase(numpy.transpose(self.recbase))
        return rec

    def cartesian(self, u):
        """return cartesian coordinates of a lattice vector"""
        rc = numpy.dot(u, self.base)
        return rc

    def dot(self, u, v):
        """return dot product of 2 lattice vectors"""
        dp = numpy.dot(u, numpy.dot(self.metrics, v))
        return dp

    def norm(self, u):
        """return norm of a lattice vector"""
        # CLF - duplicated code from dot for the sake of speed
        return math.sqrt(numpy.dot(u, numpy.dot(self.metrics, u)))

    def dist(self, u, v):
        """return distance of 2 points in lattice coordinates"""
        duv = numpy.array(u) - numpy.array(v)
        return self.norm(duv)

    def angle(self, u, v):
        """return angle(u, v) --> angle of 2 lattice vectors in degrees"""
        ca = self.dot(u, v)/( self.norm(u)*self.norm(v) )
        return math.degrees(math.acos(ca))

    def __repr__(self):
        """string representation of this lattice"""
        epsilon = 1.0e-8
        I3 = numpy.identity(3, dtype=float)
        abcABG = numpy.array([self.a, self.b, self.c,
                            self.alpha, self.beta, self.gamma] )
        rotbaseI3diff = max(numpy.reshape(numpy.fabs(self.baserot-I3), 9))
        cartlatpar = numpy.array([1.0, 1.0, 1.0 , 90.0, 90.0, 90.0])
        latpardiff = cartlatpar - self.abcABG()
        if rotbaseI3diff > epsilon:
            s = "Lattice(base=%r)" % self.base
        elif numpy.fabs(latpardiff).max() < epsilon :
            s = "Lattice()"
        else:
            s = "Lattice(a=%g, b=%g, c=%g, alpha=%g, beta=%g, gamma=%g)" % \
                    self.abcABG()
        return s

# End of Lattice

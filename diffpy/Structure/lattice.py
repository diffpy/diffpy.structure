########################################################################
#
# Structure         by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2008 trustees of the Michigan State University.
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

Module variables:
    cartesian   -- constant instance of Lattice, default Cartesian system
"""

__id__ = "$Id$"

import copy
import math
import types
import numpy
import numpy.linalg as numalg
from diffpy.Structure import LatticeError

##############################################################################
# helper functions

# cache of exact values of cosd
_EXACT_COSD = {
        0.0 : +1.0,   60.0 : +0.5,   90.0 : 0.0,  120.0 : -0.5,
      180.0 : -1.0,  240.0 : -0.5,  270.0 : 0.0,  300.0 : +0.5
}

def cosd(x):
    """Return the cosine of x (measured in degrees).
    Avoid round-off errors for exact cosine values.
    """
    return _EXACT_COSD.get(x % 360.0, math.cos(math.radians(x)))

def sind(x):
    """Return the sine of x (measured in degrees).
    Avoid round-off errors for exact sine values.
    """
    return cosd(90.0 - x)

# End of helper functions


##############################################################################
class Lattice:
    """Lattice --> general coordinate system

    Data members:

        a, b, c, alpha, beta, gamma -- read-only lattice parameters,
               unit cell angles are in degrees.  The values of lattice
               parameters are set by setLatPar() and setLatBase().
        ar, br, cr, alphar, betar, gammar -- read-only parameters of
               reciprocal lattice, set by setLatPar() and setLatBase().
               All reciprocal angles are in degrees.
        ca, cb, cg, sa, sb, sg -- read-only cosines and sines of direct
               lattice angles, they get set by setLatPar() and setLatBase()
        car, cbr, cgr, sar, sbr, sgr -- read-only cosines and sines of
               reciprocal lattice angles, set by setLatPar() and setLatBase()
        metrics  -- metrics tensor
        base     -- matrix of row base vectors in cartesian coordinates,
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
                             from specified lattice parameters.  Unit cell
                             angles are all in degrees.
        Lattice(base=abc) -- create coordinate system using the given base,
                             abc is a 3x3 matrix (or nested list), of row
                             base vectors
        Lattice(lat)      -- create a copy of existing Lattice lat
        """
        # initialize data members, they values will be set by setLatPar()
        self.a = self.b = self.c = None
        self.alpha = self.beta = self.gamma = None
        self.ca = self.cb = self.cg = None
        self.sa = self.sb = self.sg = None
        self.ar = self.br = self.cr = None
        self.alphar = self.betar = self.gammar = None
        self.car = self.cbr = self.cgr = None
        self.sar = self.sbr = self.sgr = None
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

        a, b, c, alpha, beta, gamma -- lattice parameters, unit cell angles
                    are in degrees.
        baserot  -- unit cell rotation, base = stdbase*baserot

        Note: parameters with value None will remain unchanged.

        Return self.
        """
        if a is not None: self.a = float(a)
        if b is not None: self.b = float(b)
        if c is not None: self.c = float(c)
        if alpha is not None: self.alpha = float(alpha)
        if beta is not None: self.beta = float(beta)
        if gamma is not None: self.gamma = float(gamma)
        if baserot is not None: self.baserot = numpy.array(baserot)
        (self.ca, self.sa) = (ca, sa) = ( cosd(self.alpha), sind(self.alpha) )
        (self.cb, self.sb) = (cb, sb) = ( cosd(self.beta),  sind(self.beta) )
        (self.cg, self.sg) = (cg, sg) = ( cosd(self.gamma), sind(self.gamma) )
        # Vunit is a volume of unit cell with a=b=c=1
        Vunit = math.sqrt(1.0 + 2.0*ca*cb*cg - ca*ca - cb*cb - cg*cg)
        # reciprocal lattice
        self.ar = sa/(self.a*Vunit)
        self.br = sb/(self.b*Vunit)
        self.cr = sg/(self.c*Vunit)
        self.car = car = (cb*cg - ca)/(sb*sg); sar = math.sqrt(1.0 - car**2)
        self.cbr = cbr = (ca*cg - cb)/(sa*sg); sbr = math.sqrt(1.0 - cbr**2)
        self.cgr = cgr = (ca*cb - cg)/(sa*sb); sgr = math.sqrt(1.0 - cgr**2)
        self.sar = numpy.sqrt(1.0 - car*car)
        self.sbr = numpy.sqrt(1.0 - cbr*cbr)
        self.sgr = numpy.sqrt(1.0 - cgr*cgr)
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
        """Set matrix of unit cell base vectors and calculate corresponding
        lattice parameters and stdbase, baserot and metrics tensors.

        Return self.
        """
        self.base = numpy.array(base)
        detbase = numalg.det(self.base)
        if abs(detbase) < 1.0e-8:
            emsg = "base vectors are degenerate"
            raise LatticeError(emsg)
        elif detbase < 0.0:
            emsg = "base is not right-handed"
            raise LatticeError(emsg)
        self.a = numpy.sqrt(numpy.dot(self.base[0,:], self.base[0,:]))
        self.b = numpy.sqrt(numpy.dot(self.base[1,:], self.base[1,:]))
        self.c = numpy.sqrt(numpy.dot(self.base[2,:], self.base[2,:]))
        ca = ca = numpy.dot(self.base[1,:], self.base[2,:]) / (self.b*self.c)
        cb = cb = numpy.dot(self.base[0,:], self.base[2,:]) / (self.a*self.c)
        cg = cg = numpy.dot(self.base[0,:], self.base[1,:]) / (self.a*self.b)
        sa = sa = numpy.sqrt(1.0 - ca**2)
        sb = sb = numpy.sqrt(1.0 - cb**2)
        sg = sg = numpy.sqrt(1.0 - cg**2)
        (self.ca, self.sa) = (ca, sa)
        (self.cb, self.sb) = (cb, sb)
        (self.cg, self.sg) = (cg, sg)
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
        (self.car, self.sar) = (car, sar)
        (self.cbr, self.sbr) = (cbr, sbr)
        (self.cgr, self.sgr) = (cgr, sgr)
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
        """Return a tuple of 6 lattice parameters.
        """
        rv = (self.a, self.b, self.c, self.alpha, self.beta, self.gamma)
        return rv

    def reciprocal(self):
        """Return the reciprocal lattice to self."""
        from copy import deepcopy
        rec = deepcopy(self)
        rec.setLatBase(numpy.transpose(self.recbase))
        return rec

    def cartesian(self, u):
        """return cartesian coordinates of a lattice vector"""
        rc = numpy.dot(u, self.base)
        return rc

    def fractional(self, rc):
        """return fractional coordinates of a cartesian vector"""
        u = numpy.dot(rc, self.recbase)
        return u

    def dot(self, u, v):
        """return dot product of 2 lattice vectors"""
        dp = numpy.dot(u, numpy.dot(self.metrics, v))
        return dp

    def norm(self, u):
        """return norm of a lattice vector"""
        # CLF - duplicated code from dot for the sake of speed
        return math.sqrt(numpy.dot(u, numpy.dot(self.metrics, u)))

    def dist(self, u, v):
        """Return distance of 2 points in lattice coordinates.
        """
        duv = numpy.array(u) - numpy.array(v)
        return self.norm(duv)

    def angle(self, u, v):
        """Return angle(u, v) --> angle of 2 lattice vectors in degrees.
        """
        ca = self.dot(u, v)/( self.norm(u)*self.norm(v) )
        return math.degrees(math.acos(ca))

    def __repr__(self):
        """String representation of this lattice.
        """
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

##############################################################################
# module variables

cartesian = Lattice()

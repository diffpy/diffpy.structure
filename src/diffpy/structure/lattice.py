#!/usr/bin/env python3
##############################################################################
#
# diffpy.structure  by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2008 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################

"""class Lattice stores properites and provides simple operations in lattice
coordinate system.

Module variables:
    cartesian   -- constant instance of Lattice, default Cartesian system
"""

import math
import numpy
import numpy.linalg as numalg
from diffpy.structure import LatticeError

# Helper Functions -----------------------------------------------------------

# exact values of cosd
_EXACT_COSD = {
        0.0 : +1.0,   60.0 : +0.5,   90.0 : 0.0,  120.0 : -0.5,
      180.0 : -1.0,  240.0 : -0.5,  270.0 : 0.0,  300.0 : +0.5
}

def cosd(x):
    """Return the cosine of x (measured in degrees).
    Avoid round-off errors for exact cosine values.
    """
    rv = _EXACT_COSD.get(x % 360.0)
    if rv is None:  rv = math.cos(math.radians(x))
    return rv


def sind(x):
    """Return the sine of x (measured in degrees).
    Avoid round-off errors for exact sine values.
    """
    return cosd(90.0 - x)

# ----------------------------------------------------------------------------

class Lattice:
    """
    General coordinate system and associated operations.

    Parameters
    ----------
    a : float or Lattice, optional
        The cell length *a*.  When present, other cell parameters
        must be also specified.  When of the `Lattice` type, create
        a duplicate Lattice.
    b : float
        The cell length *b*.
    c : float
        The cell length *c*.
    alpha : float
        The angle between the *b* and *c* axes in degrees.
    beta : float
        The angle between the *b* and *c* axes in degrees.
    gamma : float
        The angle between the *a* and *b* axes in degrees.
    baserot : array_like, optional
        The 3x3 rotation matrix of the base vectors with respect
        to their standard setting.
    base : array_like, optional
        The 3x3 array of row base vectors.  This must be the
        only argument when present.

    Attributes
    ----------
    metric : ndarray
        The metric tensor.
    base : ndarray
        The 3x3 matrix of row base vectors in Cartesian coordinates,
        which may be rotated, i.e., base = stdbase @ baserot.
    stdbase : ndarray
        The 3x3 matrix of row base vectors in standard orientation.
    baserot : ndarray
        The rotation matrix for the `base`.
    recbase : ndarray
        The inverse of the `base` matrix, where the columns give
        reciprocal vectors in Cartesian coordinates.
    normbase : ndarray
        The `base` vectors scaled by magnitudes of reciprocal cell lengths.
    recnormbase : ndarray
        The inverse matrix of `normbase`.
    isotropicunit : ndarray
        The 3x3 tensor for a unit isotropic displacement parameters in this
        coordinate system.  This is an identity matrix when this Lattice
        is orthonormal.

    Note
    ----
    The array attributes are read-only.  They get updated by changing
    some lattice parameters or by calling the `setLatPar()` or
    `setLatBase()` methods.

    Examples
    --------
    Create a Cartesian coordinate system::

    >>> Lattice()

    Create coordinate system with the cell lengths ``a``, ``b``, ``c``
    and cell angles ``alpha``, ``beta``, ``gamma`` in degrees::

    >>> Lattice(a, b, c, alpha, beta, gamma)

    Create a duplicate of an existing Lattice ``lat``::

    >>> Lattice(lat)

    Create coordinate system with the base vectors given by rows
    of the ``abc`` matrix::

    >>> Lattice(base=abc)
    """

    # round-off tolerance
    _epsilon = 1.0e-8


    def __init__(self, a=None, b=None, c=None,
                 alpha=None, beta=None, gamma=None,
                 baserot=None, base=None):
        # build a set of provided argument names for later use.
        apairs = (('a', a), ('b', b), ('c', c),
                  ('alpha', alpha), ('beta', beta), ('gamma', gamma),
                  ('baserot', baserot), ('base', base))
        argset = set(n for n, v in apairs if v is not None)
        # initialize data members, they values will be set by setLatPar()
        self._a = self._b = self._c = None
        self._alpha = self._beta = self._gamma = None
        self._ca = self._cb = self._cg = None
        self._sa = self._sb = self._sg = None
        self._ar = self._br = self._cr = None
        self._alphar = self._betar = self._gammar = None
        self._car = self._cbr = self._cgr = None
        self._sar = self._sbr = self._sgr = None
        self.baserot = numpy.identity(3)
        self.base = self.recbase = None
        self.normbase = self.recnormbase = None
        # work out argument variants
        # Lattice()
        if not argset:
            self.setLatPar(1.0, 1.0, 1.0, 90.0, 90.0, 90.0, baserot)
        # Lattice(base=abc)
        elif base is not None:
            if len(argset) > 1:
                raise ValueError("'base' must be the only argument.")
            self.setLatBase(base)
        # Lattice(lat)
        elif isinstance(a, Lattice):
            if len(argset) > 1:
                raise ValueError("Lattice object must be the only argument.")
            self.__dict__.update(a.__dict__)
        # otherwise do default Lattice(a, b, c, alpha, beta, gamma)
        else:
            abcabg = ('a', 'b', 'c', 'alpha', 'beta', 'gamma')
            if not argset.issuperset(abcabg):
                raise ValueError("Provide all 6 cell parameters.")
            self.setLatPar(a, b, c, alpha, beta, gamma, baserot=baserot)
        return


    def setLatPar(self, a=None, b=None, c=None,
            alpha=None, beta=None, gamma=None, baserot=None):
        """set lattice parameters and all related tensors

        Parameters
        ----------
        a : float, optional
            New value for the cell length ``a``. If not specified,
            the associated lattice parameter is thus left unchanged.
        b : float, optional
            New value for the cell length ``b``. If not specified,
            the associated lattice parameter is thus left unchanged.
        c : float, optional
            New value for the cell length ``c``. If not specified,
            the associated lattice parameter is thus left unchanged.
        alpha : float, optional
            New value for the cell angle ``alpha``. If not specified,
            the associated lattice parameter is thus left unchanged.
        beta : float, optional
            New value for the cell angle ``beta``. If not specified,
            the associated lattice parameter is thus left unchanged.
        gamma : float, optional
            New value for the cell angle ``gamma``. If not specified,
            the associated lattice parameter is thus left unchanged.
        baserot : ndarray, optional
            Unit cell rotation maxtrix, expect to be a 3x3 matrix.
        base : ndarray, optional
            Base vectors of lattice

        Note
        ----
        parameters with value None will remain unchanged.
        """
        if a is not None: self._a = float(a)
        if b is not None: self._b = float(b)
        if c is not None: self._c = float(c)
        if alpha is not None: self._alpha = float(alpha)
        if beta is not None: self._beta = float(beta)
        if gamma is not None: self._gamma = float(gamma)
        if baserot is not None: self.baserot = numpy.array(baserot)
        self._ca = ca = cosd(self.alpha)
        self._cb = cb = cosd(self.beta)
        self._cg = cg = cosd(self.gamma)
        self._sa = sa = sind(self.alpha)
        self._sb = sb = sind(self.beta)
        self._sg = sg = sind(self.gamma)
        # cache the unit volume value
        Vunit = self.unitvolume
        # reciprocal lattice
        self._ar = ar = sa/(self.a*Vunit)
        self._br = br = sb/(self.b*Vunit)
        self._cr = cr = sg/(self.c*Vunit)
        self._car = car = (cb*cg - ca)/(sb*sg)
        self._cbr = cbr = (ca*cg - cb)/(sa*sg)
        self._cgr = cgr = (ca*cb - cg)/(sa*sb)
        self._sar = math.sqrt(1.0 - car*car)
        self._sbr = math.sqrt(1.0 - cbr*cbr)
        self._sgr = sgr = math.sqrt(1.0 - cgr*cgr)
        self._alphar = math.degrees(math.acos(car))
        self._betar = math.degrees(math.acos(cbr))
        self._gammar = math.degrees(math.acos(cgr))
        # metric tensor
        self.metrics = numpy.array( [
                [ self.a*self.a,     self.a*self.b*cg,  self.a*self.c*cb ],
                [ self.b*self.a*cg,  self.b*self.b,     self.b*self.c*ca ],
                [ self.c*self.a*cb,  self.c*self.b*ca,  self.c*self.c    ] ],
                dtype=float )
        # standard Cartesian coordinates of lattice vectors
        self.stdbase = numpy.array( [
                [ 1.0/ar,    -cgr/sgr/ar,   cb*self.a ],
                [ 0.0,        self.b*sa,    self.b*ca ],
                [ 0.0,        0.0,          self.c    ] ],
                dtype=float )
        # Cartesian coordinates of lattice vectors
        self.base = numpy.dot(self.stdbase, self.baserot)
        self.recbase = numalg.inv(self.base)
        # bases normalized to unit reciprocal vectors
        self.normbase = self.base * [[ar], [br], [cr]]
        self.recnormbase = self.recbase / [ar, br, cr]
        self.isotropicunit = _isotropicunit(self.recnormbase)
        return


    def setLatBase(self, base):
        """Set matrix of unit cell base vectors and calculate corresponding
        lattice parameters and stdbase, baserot and metrics tensors.

        Parameters
        ----------
        base : ndarray, optional
            Base vectors of lattice, default to None.
        """
        self.base = numpy.array(base)
        detbase = numalg.det(self.base)
        if abs(detbase) < 1.0e-8:
            emsg = "base vectors are degenerate"
            raise LatticeError(emsg)
        elif detbase < 0.0:
            emsg = "base is not right-handed"
            raise LatticeError(emsg)
        self._a = a = math.sqrt(numpy.dot(self.base[0,:], self.base[0,:]))
        self._b = b = math.sqrt(numpy.dot(self.base[1,:], self.base[1,:]))
        self._c = c = math.sqrt(numpy.dot(self.base[2,:], self.base[2,:]))
        self._ca = ca = numpy.dot(self.base[1,:], self.base[2,:]) / (b*c)
        self._cb = cb = numpy.dot(self.base[0,:], self.base[2,:]) / (a*c)
        self._cg = cg = numpy.dot(self.base[0,:], self.base[1,:]) / (a*b)
        self._sa = sa = math.sqrt(1.0 - ca**2)
        self._sb = sb = math.sqrt(1.0 - cb**2)
        self._sg = sg = math.sqrt(1.0 - cg**2)
        self._alpha = math.degrees(math.acos(ca))
        self._beta = math.degrees(math.acos(cb))
        self._gamma = math.degrees(math.acos(cg))
        # cache the unit volume value
        Vunit = self.unitvolume
        # reciprocal lattice
        self._ar = ar = sa/(self.a*Vunit)
        self._br = br = sb/(self.b*Vunit)
        self._cr = cr = sg/(self.c*Vunit)
        self._car = car = (cb*cg - ca)/(sb*sg)
        self._cbr = cbr = (ca*cg - cb)/(sa*sg)
        self._cgr = cgr = (ca*cb - cg)/(sa*sb)
        self._sar = math.sqrt(1.0 - car**2)
        self._sbr = math.sqrt(1.0 - cbr**2)
        self._sgr = sgr = math.sqrt(1.0 - cgr**2)
        self._alphar = math.degrees(math.acos(car))
        self._betar = math.degrees(math.acos(cbr))
        self._gammar = math.degrees(math.acos(cgr))
        # standard orientation of lattice vectors
        self.stdbase = numpy.array([
                [ 1.0/ar,   -cgr/sgr/ar,    cb*a ],
                [ 0.0,       b*sa,          b*ca ],
                [ 0.0,       0.0,           c    ]],
                dtype=float)
        # calculate unit cell rotation matrix,  base = stdbase @ baserot
        self.baserot = numpy.dot(numalg.inv(self.stdbase), self.base)
        self.recbase = numalg.inv(self.base)
        # bases normalized to unit reciprocal vectors
        self.normbase = self.base * [[ar], [br], [cr]]
        self.recnormbase = self.recbase / [ar, br, cr]
        self.isotropicunit = _isotropicunit(self.recnormbase)
        # update metrics tensor
        self.metrics = numpy.array([
                [ a*a,     a*b*cg,  a*c*cb ],
                [ b*a*cg,  b*b,     b*c*ca ],
                [ c*a*cb,  c*b*ca,  c*c    ]],
                dtype=float)
        return


    def abcABG(self):
        """Return a tuple of 6 lattice parameters.

        Returns
        -------
        rv : tuple
            tuple of (a, b, c, alpha, beta, gamma)
        """
        rv = (self.a, self.b, self.c, self.alpha, self.beta, self.gamma)
        return rv


    def reciprocal(self):
        """Return the reciprocal lattice.

        Returns
        -------
        rv : Lattice
            Reciprocal lattice of current lattice.
        """
        rv = Lattice(base=numpy.transpose(self.recbase))
        return rv


    def cartesian(self, u):
        """Return Cartesian coordinates of a lattice vector.

        Parameters
        ----------
        u : ndarray
            Lattice vector in fractional coordinates or an N-by-3 array of
            lattice vectors

        Returns
        -------
        rc : ndarray
            Cartesian coordinates of a lattice vector.
        """
        rc = numpy.dot(u, self.base)
        return rc


    def fractional(self, rc):
        """Return fractional coordinates of a Cartesian vector.

        Parameters
        ----------
        rc : ndarray
            Lattice vector in Cartesian coordinates or an N-by-3 array of
            Cartesian vectors

        Returns
        -------
        u : ndarray
            fractional coordinates of a lattice vector.
        """
        u = numpy.dot(rc, self.recbase)
        return u


    def dot(self, u, v):
        """Return dot product of 2 lattice vectors.

        Parameters
        ----------
        u : ndarray
            lattice vector
        v : ndarray
            lattice vector

        Returns
        -------
        dp : dot product of two lattice vectors, u and v
        """
        dp = (u * numpy.dot(v, self.metrics)).sum(axis=-1)
        return dp


    def norm(self, xyz):
        """Calculate norm of a lattice vector.

        Parameters
        ----------
        xyz: ndarray
            vector or an N-by-3 array of fractional coordinates.

        Returns
        -------
        float or an array of the same length as xyz.
        """
        # this is a few percent faster than sqrt(dot(u, u)).
        return numpy.sqrt((self.cartesian(xyz)**2).sum(axis=-1))


    def rnorm(self, hkl):
        """Calculate norm of a reciprocal vector.

        Parameters
        ----------
        hkl : ndarray
            vector or an N-by-3 array of reciprocal coordinates.

        Returns
        -------
        float or an array of the same length as hkl.
        """
        hklcartn = numpy.dot(hkl, self.recbase.T)
        return numpy.sqrt((hklcartn**2).sum(axis=-1))


    def dist(self, u, v):
        """Calculate distance between 2 points in lattice coordinates.

        Parameters
        ----------
        u : ndarray
            vector or N-by-3 matrix of fractional coordinates.
        v : ndarray
            vector or N-by-3 matrix of fractional coordinates.

        Note
        ----
        Sizes of u, v must match when both of them are matrices.

        Returns
        -------
        float or an array of the same length as the matrix.
        """
        duv = numpy.asarray(u) - v
        return self.norm(duv)


    def angle(self, u, v):
        """Calculate angle between 2 lattice vectors in degrees.

        Parameters
        ----------
        u : ndarray
            lattice vector
        v : ndarray
            lattice vector

        Returns
        -------
        rv : float
            angle between 2 lattice vectors, u and v,  in degrees
        """
        ca = self.dot(u, v)/( self.norm(u)*self.norm(v) )
        # avoid round-off errors that would make abs(ca) greater than 1
        if numpy.isscalar(ca):
            ca = max(min(ca, 1), -1)
            rv = math.degrees(math.acos(ca))
        else:
            ca[ca < -1] = -1
            ca[ca > +1] = +1
            rv = numpy.degrees(numpy.arccos(ca))
        return rv


    def isanisotropic(self, umx):
        """Check if the atomic displacement parameter(APD) matrix
        is anisotropic in this coordinate system.

        Parameters
        ----------
        umx : ndarray
            atomic displacement parameter(adp) matrix.

        Returns
        -------
        rv : bool
            flag of  ADP matrix in this coordinate system
        """
        umx = numpy.asarray(umx)
        utr = numpy.trace(umx) / umx.shape[0]
        udmax = numpy.fabs(umx - utr * self.isotropicunit).max()
        rv = udmax > self._epsilon
        return rv


    def __repr__(self):
        """String representation of this lattice.
        """
        I3 = numpy.identity(3, dtype=float)
        rotbaseI3diff = max(numpy.reshape(numpy.fabs(self.baserot-I3), 9))
        cartlatpar = numpy.array([1.0, 1.0, 1.0 , 90.0, 90.0, 90.0])
        latpardiff = cartlatpar - self.abcABG()
        if rotbaseI3diff > self._epsilon:
            s = "Lattice(base=%r)" % self.base
        elif numpy.fabs(latpardiff).max() < self._epsilon:
            s = "Lattice()"
        else:
            s = "Lattice(a=%g, b=%g, c=%g, alpha=%g, beta=%g, gamma=%g)" % \
                    self.abcABG()
        return s

    # read-write properties

    a = property(lambda self: self._a,
                 lambda self, value: self.setLatPar(a=value),
                 doc='float: Unit cell length a')

    b = property(lambda self: self._b,
                 lambda self, value: self.setLatPar(b=value),
                 doc='float: Unit cell length b')

    c = property(lambda self: self._c,
                 lambda self, value: self.setLatPar(c=value),
                 doc='float: Unit cell length c')

    alpha = property(lambda self: self._alpha,
                     lambda self, value: self.setLatPar(alpha=value),
                     doc='float: Cell angle alpha in degrees')

    beta = property(lambda self: self._beta,
                    lambda self, value: self.setLatPar(beta=value),
                    doc='float: Cell angle beta in degrees')

    gamma = property(lambda self: self._gamma,
                     lambda self, value: self.setLatPar(gamma=value),
                     doc='float: Cell angle gamma in degrees')

    # read-only derived properties

    @property
    def unitvolume(self):
        '''Cell volume, assuming a=b=c=1
        '''
        # Recalculate lattice cosines to ensure this is right
        # even if ca, cb, cg data were not yet updated.
        ca = cosd(self.alpha)
        cb = cosd(self.beta)
        cg = cosd(self.gamma)
        rv = math.sqrt( 1.0 + 2.0*ca*cb*cg - ca*ca - cb*cb - cg*cg)
        return rv

    volume = property(lambda self: self.a * self.b * self.c * self.unitvolume,
                      doc='float: Lattice cell volume')

    ca = property(lambda self: self._ca,
                  doc='float: Cosine of the cell angle alpha')

    cb = property(lambda self: self._cb,
                  doc='float: Cosine of the cell angle beta')

    cg = property(lambda self: self._cg,
                  doc='float: Cosine of the cell angle gamma')

    sa = property(lambda self: self._sa,
                  doc='float: Sine of the cell angle alpha')

    sb = property(lambda self: self._sb,
                  doc='float: Sine of the cell angle beta')

    sg = property(lambda self: self._sg,
                  doc='float: Sine of the cell angle gamma')

    ar = property(lambda self: self._ar,
                  doc='float: Cell length a of the reciprocal lattice')

    br = property(lambda self: self._br,
                  doc='float: Cell length b of the reciprocal lattice')

    cr = property(lambda self: self._cr,
                  doc='float: Cell length c of the reciprocal lattice')

    alphar = property(lambda self: self._alphar,
                      doc='float: Reciprocal lattice angle alpha in degrees')

    betar = property(lambda self: self._betar,
                     doc='float: Reciprocal lattice angle beta in degrees')

    gammar = property(lambda self: self._gammar,
                      doc='float: Reciprocal lattice angle gamma in degrees')

    car = property(lambda self: self._car,
                   doc='float: Cosine of the reciprocal angle alpha')

    cbr = property(lambda self: self._cbr,
                   doc='float: Cosine of the reciprocal angle beta')

    cgr = property(lambda self: self._cgr,
                   doc='float: Cosine of the reciprocal angle gamma')

    sar = property(lambda self: self._sar,
                   doc='float: Sine of the reciprocal angle alpha')

    sbr = property(lambda self: self._sbr,
                   doc='flot: Sine of the reciprocal angle beta')

    sgr = property(lambda self: self._sgr,
                   doc='float: Sine of the reciprocal angle gamma')

# End of class Lattice

# Local Helpers --------------------------------------------------------------

def _isotropicunit(recnormbase):
    """Calculate matrix for unit isotropic displacement parameters.

    Parameters
    ----------
    recnormbase : ndarray
        inverse of normalized base vectors of the lattice.

    Returns
    -------
    isounit: ndarray
        matrix for unit isotropic displacement parameters
    """
    isounit = numpy.dot(recnormbase.T, recnormbase)
    # ensure there are no round-off deviations on the diagonal
    isounit[0, 0] = 1
    isounit[1, 1] = 1
    isounit[2, 2] = 1
    return isounit

# Module Constants -----------------------------------------------------------

cartesian = Lattice()

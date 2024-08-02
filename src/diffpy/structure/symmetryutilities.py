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

"""Symmetry utility functions such as expansion of asymmetric unit,
and generation of positional constraints.

Attributes
----------
epsilon : float
    Default tolerance for equality of 2 positions, also
    used for identification of special positions.

stdUsymbols : list
    Standard symbols denoting elements of anisotropic thermal
    displacement tensor.
"""

import re
import sys

import numpy

from diffpy.structure.structureerrors import SymmetryError

# Constants ------------------------------------------------------------------

epsilon = 1.0e-5

stdUsymbols = ["U11", "U22", "U33", "U12", "U13", "U23"]

# ----------------------------------------------------------------------------


def isSpaceGroupLatPar(spacegroup, a, b, c, alpha, beta, gamma):
    """Check if space group allows passed lattice parameters.

    Parameters
    ----------
    spacegroup : SpaceGroup
        Instance of `SpaceGroup`.
    a, b, c, alpha, beta, gamma : float
        `Lattice` parameters.

    Return
    ------
    bool
        ``True`` when lattice parameters are allowed by space group.

    Note
    ----
    Crystal system rules:

    Benjamin, W. A., Introduction to crystallography, New York (1969), p.60.
    """

    # crystal system rules
    # ref: Benjamin, W. A., Introduction to crystallography,
    # New York (1969), p.60
    def check_triclinic():
        return True

    def check_monoclinic():
        rv = (alpha == gamma == 90) or (alpha == beta == 90)
        return rv

    def check_orthorhombic():
        return alpha == beta == gamma == 90

    def check_tetragonal():
        return a == b and alpha == beta == gamma == 90

    def check_trigonal():
        rv = (a == b == c and alpha == beta == gamma) or (a == b and alpha == beta == 90 and gamma == 120)
        return rv

    def check_hexagonal():
        return a == b and alpha == beta == 90 and gamma == 120

    def check_cubic():
        return a == b == c and alpha == beta == gamma == 90

    crystal_system_rules = {
        "TRICLINIC": check_triclinic,
        "MONOCLINIC": check_monoclinic,
        "ORTHORHOMBIC": check_orthorhombic,
        "TETRAGONAL": check_tetragonal,
        "TRIGONAL": check_trigonal,
        "HEXAGONAL": check_hexagonal,
        "CUBIC": check_cubic,
    }
    rule = crystal_system_rules[spacegroup.crystal_system]
    return rule()


# Constant regular expression used in isconstantFormula().
# isconstantFormula runs faster when regular expression is not
# compiled per every single call.

_rx_constant_formula = re.compile(r"[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)??(/[-+]?\d+)?$")


def isconstantFormula(s):
    """Check if formula string is constant.

    Parameters
    ----------
    s : str
        Formula string.

    Return
    ------
    bool
        ``True`` when argument is a floating point number or a fraction of float with integer.
    """
    res = _rx_constant_formula.match(s.replace(" ", ""))
    return bool(res)


# Helper class intended for this module only:


class _Position2Tuple(object):
    """Create callable object that converts fractional coordinates to
    a tuple of integers with given precision. For presision close to zero
    it will return a tuples of double.

    Note
    ----
    Helper class intended for local use only.

    Parameters
    ----------
    eps : float
        Cutoff for equivalent coordinates.

    Attributes
    ----------
    eps : float
        Cutoff for equivalent coordinates. When two coordiantes map to the
        same tuple, they are closer than `eps`.
    """

    def __init__(self, eps=None):
        if eps is None:
            eps = epsilon
        # ensure self.eps has exact machine representation
        self.eps = eps + 1.0
        self.eps = self.eps - 1.0
        # no conversions for very small eps
        if self.eps == 0.0 or 1.0 / self.eps > sys.maxsize:
            self.eps = 0.0
        return

    def __call__(self, xyz):
        """Convert array of fractional coordinates to a tuple.

        Parameters
        ----------
        xyz : Iterable
            Fractional coordinates.

        Return
        ------
        tuple
            Tuple of 3 float when `eps` is zero, otherwise tuple of 3 int.
        """
        # no conversion case
        if self.eps == 0.0:
            tpl = tuple(xyz % 1.0)
            return tpl
        # here we convert to integer
        tpl = tuple([int((xi - numpy.floor(xi)) / self.eps) for xi in xyz])
        return tpl


# End of class _Position2Tuple


def positionDifference(xyz0, xyz1):
    """Smallest difference between two coordinates in periodic lattice.

    Parameters
    ----------
    xyz0, xyz1 : array_like
        Fractional coordinates.

    Return
    ------
    dxyz : numpy.ndarray
        Smallest difference between two coordinates in periodic lattice
        with ``0 <= dxyz <= 0.5``.
    """
    dxyz = numpy.asarray(xyz0) - xyz1
    # map differences to [0,0.5]
    dxyz = dxyz - numpy.floor(dxyz)
    mask = dxyz > 0.5
    dxyz[mask] = 1.0 - dxyz[mask]
    return dxyz


def nearestSiteIndex(sites, xyz):
    """Index of the nearest site to a specified position.

    Parameters
    ----------
    sites : array_like
        List of coordinates.
    xyz : array_like
        Single position.

    Return
    ------
    int
        Index of the nearest site.
    """
    # we use box distance to be consistent with _Position2Tuple conversion
    dbox = positionDifference(sites, xyz).max(axis=1)
    nearindex = numpy.argmin(dbox)
    return nearindex


def equalPositions(xyz0, xyz1, eps):
    """Equality of two coordinates with optional tolerance.

    Parameters
    ----------
    xyz0, xyz1 : array_like
        Fractional coordinates.
    eps : float
        Tolerance for equality of coordinates.

    Return
    ------
    bool
        ``True`` when two coordinates are closer than `eps`.
    """
    # we use box distance to be consistent with _Position2Tuple conversion
    dxyz = positionDifference(xyz0, xyz1)
    return numpy.all(dxyz <= eps)


def expandPosition(spacegroup, xyz, sgoffset=[0, 0, 0], eps=None):
    """Obtain unique equivalent positions and corresponding operations.

    Parameters
    ----------
    spacegroup : SpaceGroup
        Instance of SpaceGroup.
    xyz : list
        Position to be expanded.
    sgoffset : list, Optional
        Offset of space group origin ``[0, 0, 0]``. Default is ``[0, 0, 0]``.
    eps : float, Optional
        Cutoff for equal positions, default is ``1.0e-5``.

    Return
    ------
    tuple
        A tuple with ``(list of unique equivalent positions, nested
        list of `SpaceGroups.SymOp` instances, site multiplicity)``.
    """
    sgoffset = numpy.asarray(sgoffset, dtype=float)
    if eps is None:
        eps = epsilon
    pos2tuple = _Position2Tuple(eps)
    positions = []
    site_symops = {}  # position tuples with [related symops]
    for symop in spacegroup.iter_symops():
        # operate on coordinates in non-shifted spacegroup
        pos = symop(xyz + sgoffset) - sgoffset
        mask = numpy.logical_or(pos < 0.0, pos >= 1.0)
        pos[mask] -= numpy.floor(pos[mask])
        tpl = pos2tuple(pos)
        if tpl not in site_symops:
            pos_is_new = True
            site_symops[tpl] = []
            # double check if there is any position nearby
            if positions:
                nearpos = positions[nearestSiteIndex(positions, pos)]
                # is it an equivalent position?
                if equalPositions(nearpos, pos, eps):
                    # tpl should map to the same list as nearpos
                    site_symops[tpl] = site_symops[pos2tuple(nearpos)]
                    pos_is_new = False
            if pos_is_new:
                positions.append(pos)
        # here tpl is inside site_symops
        site_symops[tpl].append(symop)
    # pos_symops is nested list of symops associated with each position
    pos_symops = [site_symops[pos2tuple(p)] for p in positions]
    multiplicity = len(positions)
    return positions, pos_symops, multiplicity


def nullSpace(A):
    """Null space of matrix A."""
    from numpy import linalg

    u, s, v = linalg.svd(A)
    # s may have smaller dimension than v
    vnrows = numpy.shape(v)[0]
    mask = numpy.ones(vnrows, dtype=bool)
    mask[s > epsilon] = False
    null_space = numpy.compress(mask, v, axis=0)
    return null_space


def _findInvariants(symops):
    """Find a list of symmetry operations which contains identity.

    Parameters
    ----------
    symops : list of SymOp
        Nested list of `SymOp` instances.

    Return
    ------
    list
        List-item in symops which contains identity.

    Raise
    -----
    ValueError
        When identity is not found.
    """
    invrnts = None
    R0 = numpy.identity(3, dtype=float)
    t0 = numpy.zeros(3, dtype=float)
    for ops in symops:
        for op in ops:
            if numpy.all(op.R == R0) and numpy.all(op.t == t0):
                invrnts = ops
                break
        if invrnts:
            break
    if invrnts is None:
        emsg = "Could not find identity operation."
        raise ValueError(emsg)
    return invrnts


# ----------------------------------------------------------------------------


class GeneratorSite(object):
    """Storage of data related to a generator positions.

    Parameters
    ----------
    spacegroup : SpaceGroup
        Instance of `SpaceGroup`.
    xyz : array_like
        Generating site. When `xyz` is close to special
        position `self.xyz` will be adjusted.
    Uij : array_like, Optional
        Thermal factors at generator site. Yields `self.Uij`
        after adjusting to spacegroup symmetry. Default is zeros.
    sgoffset : list, Optional
        Offset of space group origin ``[0, 0, 0]``. Default is ``[0, 0, 0]``.
    eps : float, Optional
        Cutoff for equal positions. Default is ``1.0e-5``.

    Attributes
    ----------
    xyz : numpy.ndarray
        Fractional coordinates of generator site.
    Uij : numpy.ndarray
        Anisotropic thermal displacement at generator site.
    sgoffset : numpy.ndarray
        Offset of space group origin ``[0, 0, 0]``.
    eps : float
        Cutoff for equal positions.
    eqxyz : list
        List of equivalent positions.
    eqUij : list
        List of displacement matrices at equivalent positions.
    symops : list
        Nested list of operations per each `eqxyz`.
    multiplicity : int
        Generator site multiplicity.
    Uisotropy : bool
        Bool flag for isotropic thermal factors.
    invariants : list
        List of invariant operations for generator site.
    null_space : numpy.ndarray
        Null space of all possible differences of invariant
        rotational matrices, this is a base of symmetry
        allowed shifts.
    Uspace : numpy.ndarray
        3D array of independent components of U matrices.
    pparameters : list
        List of ``(xyz symbol, value)`` pairs.
    Uparameters : list
        List of ``(U symbol, value)`` pairs.
    """

    Ucomponents = numpy.array(
        [
            [[1, 0, 0], [0, 0, 0], [0, 0, 0]],
            [[0, 0, 0], [0, 1, 0], [0, 0, 0]],
            [[0, 0, 0], [0, 0, 0], [0, 0, 1]],
            [[0, 1, 0], [1, 0, 0], [0, 0, 0]],
            [[0, 0, 1], [0, 0, 0], [1, 0, 0]],
            [[0, 0, 0], [0, 0, 1], [0, 1, 0]],
        ],
        dtype=float,
    )
    """numpy.ndarray: 6x3x3 array of independent components of U matrices."""

    idx2Usymbol = {0: "U11", 1: "U12", 2: "U13", 3: "U12", 4: "U22", 5: "U23", 6: "U13", 7: "U23", 8: "U33"}
    """dict: Mapping of index to standard U symbol."""

    def __init__(self, spacegroup, xyz, Uij=numpy.zeros((3, 3)), sgoffset=[0, 0, 0], eps=None):
        if eps is None:
            eps = epsilon
        # just declare the members
        self.xyz = numpy.array(xyz, dtype=float)
        self.Uij = numpy.array(Uij, dtype=float)
        self.sgoffset = numpy.array(sgoffset, dtype=float)
        self.eps = eps
        self.eqxyz = []
        self.eqUij = []
        self.symops = None
        self.multiplicity = None
        self.Uisotropy = False
        self.invariants = []
        self.null_space = None
        self.Uspace = None
        self.pparameters = []
        self.Uparameters = []
        # fill in the values
        sites, ops, mult = expandPosition(spacegroup, xyz, sgoffset, eps)
        invariants = _findInvariants(ops)
        # shift self.xyz exactly to the special position
        if mult > 1:
            xyzdups = numpy.array([op(xyz + self.sgoffset) - self.sgoffset for op in invariants])
            dxyz = xyzdups - xyz
            dxyz = numpy.mean(dxyz - dxyz.round(), axis=0)
            # recalculate if needed
            if numpy.any(dxyz != 0.0):
                self.xyz = xyz + dxyz
                self.xyz[numpy.fabs(self.xyz) < self.eps] = 0.0
                sites, ops, mult = expandPosition(spacegroup, self.xyz, self.sgoffset, eps)
                invariants = _findInvariants(ops)
        # self.xyz, sites, ops are all adjusted here
        self.eqxyz = sites
        self.symops = ops
        self.multiplicity = mult
        self.invariants = invariants
        self._findNullSpace()
        self._findPosParameters()
        self._findUSpace()
        self._findUParameters()
        self._findeqUij()
        return

    def signedRatStr(self, x):
        """Convert floating point number to signed rational representation.

        Possible fractional are multiples of 1/3, 1/6, 1/7, 1/9, if these
        are not close, return `%+g` format.

        Parameters
        ----------
        x : float
            Floating point number.

        Return
        ------
        str
            Signed rational representation of `x`.
        """
        s = "{:.8g}".format(x)
        if len(s) < 6:
            return "%+g" % x
        den = numpy.array([3.0, 6.0, 7.0, 9.0])
        nom = x * den
        idx = numpy.where(numpy.fabs(nom - nom.round()) < self.eps)[0]
        if idx.size == 0:
            return "%+g" % x
        # here we have fraction
        return "%+.0f/%.0f" % (nom[idx[0]], den[idx[0]])

    def _findNullSpace(self):
        """Calculate `self.null_space` from `self.invariants`.
        Try to represent `self.null_space` using small integers.
        """
        R0 = self.invariants[0].R
        Rdiff = [(symop.R - R0) for symop in self.invariants]
        Rdiff = numpy.concatenate(Rdiff, axis=0)
        self.null_space = nullSpace(Rdiff)
        if self.null_space.size == 0:
            return
        # reverse sort rows of null_space rows by absolute value
        key = tuple(numpy.fabs(numpy.transpose(self.null_space))[::-1])
        order = numpy.lexsort(key)
        self.null_space = self.null_space[order[::-1]]
        # rationalize by the smallest element larger than cutoff
        cutoff = 1.0 / 32
        for row in self.null_space:
            abrow = numpy.abs(row)
            sgrow = numpy.sign(row)
            # equalize items with round-off-equal absolute value
            ii = abrow.argsort()
            delta = 1e-8 * abrow[ii[-1]]
            for k in ii[1:]:
                if abrow[k] - abrow[k - 1] < delta:
                    abrow[k] = abrow[k - 1]
            # find the smallest nonzero absolute element
            jnz = numpy.flatnonzero(abrow > cutoff)
            idx = jnz[abrow[jnz].argmin()]
            row[:] = (sgrow * abrow) / sgrow[idx] / abrow[idx]
        return

    def _findPosParameters(self):
        """Find pparameters and their values for expressing `self.xyz`."""
        usedsymbol = {}
        # parameter values depend on offset of self.xyz
        txyz = self.xyz
        # define txyz such that most of its elements are zero
        for nvec in self.null_space:
            idx = numpy.where(numpy.fabs(nvec) >= epsilon)[0][0]
            varvalue = txyz[idx] / nvec[idx]
            txyz = txyz - varvalue * nvec
            # determine standard parameter name
            vname = [s for s in "xyz"[idx:] if s not in usedsymbol][0]
            self.pparameters.append((vname, varvalue))
            usedsymbol[vname] = True
        return

    def _findUSpace(self):
        """Find independent U components with respect to invariant
        rotations.
        """
        n = len(self.invariants)
        R6zall = numpy.tile(-numpy.identity(6, dtype=float), (n, 1))
        R6zall_iter = numpy.split(R6zall, n, axis=0)
        i6kl = ((0, (0, 0)), (1, (1, 1)), (2, (2, 2)), (3, (0, 1)), (4, (0, 2)), (5, (1, 2)))
        for op, R6z in zip(self.invariants, R6zall_iter):
            R = op.R
            for j, Ucj in enumerate(self.Ucomponents):
                Ucj2 = numpy.dot(R, numpy.dot(Ucj, R.T))
                for i, kl in i6kl:
                    R6z[i, j] += Ucj2[kl]
        Usp6 = nullSpace(R6zall)
        # normalize Usp6 by its maximum component
        mxcols = numpy.argmax(numpy.fabs(Usp6), axis=1)
        mxrows = numpy.arange(len(mxcols))
        Usp6 /= Usp6[mxrows, mxcols].reshape(-1, 1)
        Usp6 = numpy.around(Usp6, 2)
        # normalize again after rounding to get correct signs
        mxcols = numpy.argmax(numpy.fabs(Usp6), axis=1)
        Usp6 /= Usp6[mxrows, mxcols].reshape(-1, 1)
        self.Uspace = numpy.tensordot(Usp6, self.Ucomponents, axes=(1, 0))
        self.Uisotropy = len(self.Uspace) == 1
        return

    def _findUParameters(self):
        """Find Uparameters and their values for expressing `self.Uij`."""
        # permute indices as     00 11 22 01 02 12 10 20 21
        diagorder = numpy.array((0, 4, 8, 1, 2, 5, 3, 6, 7))
        Uijflat = self.Uij.flatten()
        for Usp in self.Uspace:
            Uspflat = Usp.flatten()
            Uspnorm2 = numpy.dot(Uspflat, Uspflat)
            permidx = next(i for i, x in enumerate(Uspflat[diagorder]) if x == 1)
            idx = diagorder[permidx]
            vname = self.idx2Usymbol[idx]
            varvalue = numpy.dot(Uijflat, Uspflat) / Uspnorm2
            self.Uparameters.append((vname, varvalue))
        return

    def _findeqUij(self):
        """Adjust `self.Uij` and `self.eqUij` to be consistent with spacegroup."""
        self.Uij = numpy.zeros((3, 3), dtype=float)
        for i in range(len(self.Uparameters)):
            Usp = self.Uspace[i]
            varvalue = self.Uparameters[i][1]
            self.Uij += varvalue * Usp
        # now determine eqUij
        for ops in self.symops:
            # take first rotation matrix
            R = ops[0].R
            Rt = R.transpose()
            self.eqUij.append(numpy.dot(R, numpy.dot(self.Uij, Rt)))
        return

    def positionFormula(self, pos, xyzsymbols=("x", "y", "z")):
        """Formula of equivalent position with respect to generator site.

        Parameters
        ----------
        pos : array_like
            Fractional coordinates of possibly equivalent site.
        xyzsymbols : tuple, Optional
            Symbols for parametrized coordinates.

        Return
        ------
        dict
            Position formulas in a dictionary with keys equal ``("x", "y", "z")``
            or an empty dictionary when pos is not equivalent to generator.
            Formulas are formatted as ``[[-][%g*]{x|y|z}] [{+|-}%g]``, for example
            ``-x``, ``z +0.5``, ``0.25``.
        """
        # find pos in eqxyz
        idx = nearestSiteIndex(self.eqxyz, pos)
        eqpos = self.eqxyz[idx]
        if not equalPositions(eqpos, pos, self.eps):
            return {}
        # any rotation matrix should do fine
        R = self.symops[idx][0].R
        nsrotated = numpy.dot(self.null_space, numpy.transpose(R))
        # build formulas using eqpos
        # find offset
        teqpos = numpy.array(eqpos)
        for nvec, (vname, varvalue) in zip(nsrotated, self.pparameters):
            teqpos -= nvec * varvalue
        # map varnames to xyzsymbols
        name2sym = dict(zip(("x", "y", "z"), xyzsymbols))
        xyzformula = 3 * [""]
        for nvec, (vname, ignore) in zip(nsrotated, self.pparameters):
            for i in range(3):
                if abs(nvec[i]) < epsilon:
                    continue
                xyzformula[i] += "%s*%s " % (self.signedRatStr(nvec[i]), name2sym[vname])
        # add constant offset teqpos to all formulas
        for i in range(3):
            if xyzformula[i] and abs(teqpos[i]) < epsilon:
                continue
            xyzformula[i] += self.signedRatStr(teqpos[i])
        # reduce unnecessary +1* and -1*
        xyzformula = [re.sub("^[+]1[*]|(?<=[+-])1[*]", "", f).strip() for f in xyzformula]
        return dict(zip(("x", "y", "z"), xyzformula))

    def UFormula(self, pos, Usymbols=stdUsymbols):
        """List of atom displacement formulas with custom parameter symbols.

        Parameters
        ----------
        pos : array_like
            Fractional coordinates of possibly equivalent site.
        Usymbols : list, Optional
            6 symbols for possible U matrix parameters, default is
            ``["U11", "U22", "U33", "U12", "U13", "U23"]``.

        Return
        ------
        Uformula : dict
            U element formulas in a dictionary where keys are from
            ``('U11','U22','U33','U12','U13','U23')`` or empty dictionary when
            pos is not equivalent to generator.
        """
        # find pos in eqxyz
        idx = nearestSiteIndex(self.eqxyz, pos)
        eqpos = self.eqxyz[idx]
        if not equalPositions(eqpos, pos, self.eps):
            return {}
        # any rotation matrix should do fine
        R = self.symops[idx][0].R
        Rt = R.transpose()
        Usrotated = [numpy.dot(R, numpy.dot(Us, Rt)) for Us in self.Uspace]
        Uformula = dict.fromkeys(stdUsymbols, "")
        name2sym = dict(zip(stdUsymbols, Usymbols))
        for Usr, (vname, ignore) in zip(Usrotated, self.Uparameters):
            # avoid adding off-diagonal elements twice
            assert numpy.all(Usr == Usr.T)
            Usr -= numpy.tril(Usr, -1)
            Usrflat = Usr.flatten()
            for i in numpy.where(Usrflat)[0]:
                f = "%+g*%s" % (Usrflat[i], name2sym[vname])
                smbl = self.idx2Usymbol[i]
                Uformula[smbl] += f
        for smbl, f in Uformula.items():
            if not f:
                f = "0"
            f = re.sub(r"^[+]?1[*]|^[+](?=\d)|(?<=[+-])1[*]", "", f).strip()
            Uformula[smbl] = f
        return Uformula

    def eqIndex(self, pos):
        """Index of the nearest generator equivalent site.

        Parameters
        ----------
        pos : array_like
            Fractional coordinates.

        Return
        ------
        int
            Index of the nearest generator equivalent site.
        """
        return nearestSiteIndex(self.eqxyz, pos)


# End of class GeneratorSite

# ----------------------------------------------------------------------------


class ExpandAsymmetricUnit(object):
    """Expand asymmetric unit and anisotropic thermal displacement.

    Parameters
    ----------
    spacegroup : SpaceGroup
        Instance of `SpaceGroup`.
    corepos : array_like
        List of positions in asymmetric unit,
        it may contain duplicates.
    coreUijs : numpy.ndarray, Optional
        Thermal factors for `corepos`.
    sgoffset : list, Optional
        Offset of space group origin ``[0, 0, 0]``. Default is ``[0, 0, 0]``.
    eps : float, Optional
        Cutoff for duplicate positions. Default is ``1.0e-5``.

    Attributes
    ----------
    spacegroup : SpaceGroup
        Instance of `SpaceGroup`.
    corepos : array_like
        List of positions in asymmetric unit,
        it may contain duplicates.
    coreUijs : numpy.ndarray
        Thermal factors for `corepos`. Defaults to zeros.
    sgoffset : numpy.ndarray
        Offset of space group origin ``[0, 0, 0]``. Default to zeros.
    eps : float
        Cutoff for equivalent positions. Default is ``1.0e-5``.
    multiplicity : list
        Multiplicity of each site in `corepos`.
    Uisotropy : list
        Bool flags for isotropic sites in `corepos`.
    expandedpos : list
        List of equivalent positions per each site in `corepos`.
    expandedUijs : list
        List of thermal factors per each site in `corepos`.
    """

    # By design Atom instances are not accepted as arguments to keep
    # number of required imports low.
    def __init__(self, spacegroup, corepos, coreUijs=None, sgoffset=[0, 0, 0], eps=None):
        if eps is None:
            eps = epsilon
        # declare data members
        self.spacegroup = spacegroup
        self.corepos = corepos
        self.coreUijs = None
        self.sgoffset = numpy.array(sgoffset)
        self.eps = eps
        self.multiplicity = []
        self.Uisotropy = []
        self.expandedpos = []
        self.expandedUijs = []
        # obtain their values
        corelen = len(self.corepos)
        if coreUijs:
            self.coreUijs = coreUijs
        else:
            self.coreUijs = numpy.zeros((corelen, 3, 3), dtype=float)
        for cpos, cUij in zip(self.corepos, self.coreUijs):
            gen = GeneratorSite(self.spacegroup, cpos, cUij, self.sgoffset, self.eps)
            self.multiplicity.append(gen.multiplicity)
            self.Uisotropy.append(gen.Uisotropy)
            self.expandedpos.append(gen.eqxyz)
            self.expandedUijs.append(gen.eqUij)
        return


# End of class ExpandAsymmetricUnit


# Helper function for SymmetryConstraints class. It may be useful
# elsewhere therefore its name does not start with underscore.


def pruneFormulaDictionary(eqdict):
    """Remove constant items from formula dictionary.

    Parameters
    ----------
    eqdict : dict
        Formula dictionary which maps standard variable symbols
        ``("x", "U11")`` to string formulas ``("0", "-x3", "z7 +0.5")``.

    Return
    ------
    dict
        Pruned formula dictionary.
    """
    pruned = {}
    for smb, eq in eqdict.items():
        if not isconstantFormula(eq):
            pruned[smb] = eq
    return pruned


class SymmetryConstraints(object):
    """Generate symmetry constraints for specified positions.

    Parameters
    ----------
    spacegroup : SpaceGroup
        Instance of `SpaceGroup`.
    positions : array_like
        List of all positions to be constrained.
    Uijs : array_like, Optional
        List of U matrices for all constrained positions.
    sgoffset : list, Optional
        Offset of space group origin ``[0, 0, 0]``. Default is ``[0, 0, 0]``.
    eps : float, Optional
        Cutoff for duplicate positions. Default is ``1.0e-5``.

    Attributes
    ----------
    spacegroup : SpaceGroup
        Instance of `SpaceGroup`.
    positions : numpy.ndarray
        All positions to be constrained.
    Uijs : numpy.ndarray
        Thermal factors for all positions. Defaults to zeros.
    sgoffset : numpy.ndarray
        Optional offset of space group origin ``[0, 0, 0]``.
    eps : float
        Cutoff for equivalent positions. Default is ``1.0e-5``.
    corepos : list
        List of of positions in the asymmetric unit.
    coremap : dict
        Dictionary mapping indices of asymmetric core positions
        to indices of all symmetry related positions.
    poseqns : list
        List of coordinate formula dictionaries per each site.
        Formula dictionary keys are from ``("x", "y", "z")`` and
        the values are formatted as ``[[-]{x|y|z}%i] [{+|-}%g]``,
        for example: ``x0``, ``-x3``, ``z7 +0.5``, ``0.25``.
    pospars : list
        List of ``(xyz symbol, value)`` pairs.
    Ueqns : list
        List of anisotropic atomic displacement formula
        dictionaries per each position. Formula dictionary
        keys are from ``('U11','U22','U33','U12','U13','U23')``
        and the values are formatted as ``{[%g*][Uij%i]|0}``,
        for example: ``U110``, ``0.5*U2213``, ``0``.
    Upars : list
        List of ``(U symbol, value)`` pairs.
    Uisotropy : list
        List of bool flags for isotropic thermal displacements.
    """

    def __init__(self, spacegroup, positions, Uijs=None, sgoffset=[0, 0, 0], eps=None):
        if eps is None:
            eps = epsilon
        # fill in data members
        self.spacegroup = spacegroup
        self.positions = None
        self.Uijs = None
        self.sgoffset = numpy.array(sgoffset)
        self.eps = eps
        self.corepos = []
        self.coremap = {}
        self.poseqns = None
        self.pospars = []
        self.Ueqns = None
        self.Upars = []
        self.Uisotropy = None
        # handle list of lists returned by ExpandAsymmetricUnit
        if len(positions) and isinstance(positions[0], list):
            # concatenate lists before converting to Nx3 array
            flatpos = sum(positions, [])
            flatpos = numpy.array(flatpos, dtype=float).flatten()
            self.positions = flatpos.reshape((-1, 3))
        # otherwise convert to array
        else:
            flatpos = numpy.array(positions, dtype=float).flatten()
            self.positions = flatpos.reshape((-1, 3))
        # here self.positions should be a 2D numpy array
        numpos = len(self.positions)
        # adjust Uijs if not specified
        if Uijs is not None:
            self.Uijs = numpy.array(Uijs, dtype=float)
        else:
            self.Uijs = numpy.zeros((numpos, 3, 3), dtype=float)
        self.poseqns = numpos * [None]
        self.Ueqns = numpos * [None]
        self.Uisotropy = numpos * [False]
        # all members should be initialized here
        self._findConstraints()
        return

    def _findConstraints(self):
        """Find constraints for positions and anisotropic displacements `Uij`."""
        numpos = len(self.positions)
        # canonical xyzsymbols and Usymbols
        xyzsymbols = [smbl + str(i) for i in range(numpos) for smbl in "xyz"]
        Usymbols = [smbl + str(i) for i in range(numpos) for smbl in stdUsymbols]
        independent = set(range(numpos))
        for genidx in range(numpos):
            if genidx not in independent:
                continue
            # it is a generator
            self.coremap[genidx] = []
            genpos = self.positions[genidx]
            genUij = self.Uijs[genidx]
            gen = GeneratorSite(self.spacegroup, genpos, genUij, self.sgoffset, self.eps)
            # append new pparameters if there are any
            gxyzsymbols = xyzsymbols[3 * genidx : 3 * (genidx + 1)]
            for k, v in gen.pparameters:
                smbl = gxyzsymbols["xyz".index(k)]
                self.pospars.append((smbl, v))
            gUsymbols = Usymbols[6 * genidx : 6 * (genidx + 1)]
            for k, v in gen.Uparameters:
                smbl = gUsymbols[stdUsymbols.index(k)]
                self.Upars.append((smbl, v))
            # search for equivalents inside indies
            indies = sorted(independent)
            for indidx in indies:
                indpos = self.positions[indidx]
                formula = gen.positionFormula(indpos, gxyzsymbols)
                # formula is empty when indidx is independent
                if not formula:
                    continue
                # indidx is dependent here
                independent.remove(indidx)
                self.coremap[genidx].append(indidx)
                self.poseqns[indidx] = formula
                self.Ueqns[indidx] = gen.UFormula(indpos, gUsymbols)
                # make sure positions and Uijs are consistent with spacegroup
                eqidx = gen.eqIndex(indpos)
                dxyz = gen.eqxyz[eqidx] - indpos
                self.positions[indidx] += dxyz - dxyz.round()
                self.Uijs[indidx] = gen.eqUij[eqidx]
                self.Uisotropy[indidx] = gen.Uisotropy
        # all done here
        coreidx = sorted(self.coremap.keys())
        self.corepos = [self.positions[i] for i in coreidx]
        return

    def posparSymbols(self):
        """Return list of standard position parameter symbols."""
        return [n for n, v in self.pospars]

    def posparValues(self):
        """Return list of position parameters values."""
        return [v for n, v in self.pospars]

    def positionFormulas(self, xyzsymbols=None):
        """List of position formulas with custom parameter symbols.

        Parameters
        ----------
        xyzsymbols : list, Optional
            List of custom symbols used in formula strings.

        Return
        ------
        list
            List of coordinate formulas dictionaries. Formulas dictionary
            keys are from ``("x", "y", "z")`` and the values are formatted as
            ``[[-]{symbol}] [{+|-}%g]``, for example: ``x0``, ``-sym``, ``@7 +0.5``, ``0.25``.
        """
        if not xyzsymbols:
            return list(self.poseqns)
        # check xyzsymbols
        if len(xyzsymbols) < len(self.pospars):
            emsg = "Not enough symbols for %i position parameters" % len(self.pospars)
            raise SymmetryError(emsg)
        # build translation dictionary
        trsmbl = dict(zip(self.posparSymbols(), xyzsymbols))

        def translatesymbol(matchobj):
            return trsmbl[matchobj.group(0)]

        pat = re.compile(r"\b[xyz]\d+")
        rv = []
        for eqns in self.poseqns:
            treqns = {}
            for smbl, eq in eqns.items():
                treqns[smbl] = re.sub(pat, translatesymbol, eq)
            rv.append(treqns)
        return rv

    def positionFormulasPruned(self, xyzsymbols=None):
        """List of position formula dictionaries with constant items removed.

        See also
        --------
        positionFormulas()

        Parameters
        ----------
        xyzsymbols : list, Optional
            List of custom symbols used in formula strings.

        Return
        ------
        list
            List of coordinate formula dictionaries.
        """
        rv = [pruneFormulaDictionary(eqns) for eqns in self.positionFormulas(xyzsymbols)]
        return rv

    def UparSymbols(self):
        """Return list of standard atom displacement parameter symbols."""
        return [n for n, v in self.Upars]

    def UparValues(self):
        """Return list of atom displacement parameters values."""
        return [v for n, v in self.Upars]

    def UFormulas(self, Usymbols=None):
        """List of atom displacement formulas with custom parameter symbols.

        Parameters
        ----------
        Usymbols : list, Optional
            List of custom symbols used in formula strings.

        Return
        ------
        list
            List of atom displacement formula dictionaries per each site.
            Formula dictionary keys are from ``('U11','U22','U33','U12','U13','U23')``
            and the values are formatted as ``{[%g*][Usymbol]|0}``, for example:
            ``U11``, ``0.5*@37``, ``0``.
        """
        if not Usymbols:
            return list(self.Ueqns)
        # check Usymbols
        if len(Usymbols) < len(self.Upars):
            emsg = "Not enough symbols for %i U parameters" % len(self.Upars)
            raise SymmetryError(emsg)
        # build translation dictionary
        trsmbl = dict(zip(self.UparSymbols(), Usymbols))

        def translatesymbol(matchobj):
            return trsmbl[matchobj.group(0)]

        pat = re.compile(r"\bU\d\d\d+")
        rv = []
        for eqns in self.Ueqns:
            treqns = {}
            for smbl, eq in eqns.items():
                treqns[smbl] = re.sub(pat, translatesymbol, eq)
            rv.append(treqns)
        return rv

    def UFormulasPruned(self, Usymbols=None):
        """List of atom displacement formula dictionaries with constant items
        removed.

        See Also
        --------
        UFormulas()

        Parameters
        ----------
        Usymbols : list, Optional
            List of custom symbols used in formula strings.

        Return
        ------
        list
            List of atom displacement formulas in tuples of
            ``(U11, U22, U33, U12, U13, U23)``.
        """
        rv = [pruneFormulaDictionary(eqns) for eqns in self.UFormulas(Usymbols)]
        return rv


# End of class SymmetryConstraints

# ----------------------------------------------------------------------------

# basic demonstration
if __name__ == "__main__":
    from diffpy.structure.spacegroups import sg100

    site = [0.125, 0.625, 0.13]
    Uij = [[1, 2, 3], [2, 4, 5], [3, 5, 6]]
    g = GeneratorSite(sg100, site, Uij=Uij)
    fm100 = g.positionFormula(site)
    print("g = GeneratorSite(sg100, %r)" % site)
    print("g.positionFormula(%r) = %s" % (site, fm100))
    print("g.pparameters =", g.pparameters)
    print("g.Uparameters =", g.Uparameters)
    print("g.UFormula(%r) =" % site, g.UFormula(site))

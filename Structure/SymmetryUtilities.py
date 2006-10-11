"""Symmetry utility functions such as expansion of asymmetric unit,
and generation of positional constraints.
"""

# version
__id__ = '$Id$'

import sys
import re
import numpy

# Constants:
# Default tolerance for equality of 2 positions, also
# used for identification of special positions.
epsilon = 1.0e-5

def isSpaceGroupLatPar(spacegroup, a, b, c, alpha, beta, gamma):
    """Check if space group allows passed lattice parameters

    spacegroup -- instance of SpaceGroup
    a, b, c, alpha, beta, gamma -- lattice parameters

    Return bool
    """
    # ref: Benjamin, W. A., Introduction to crystallography, 
    # New York (1969), p.60
    crystal_system_rules = {
      "TRICLINIC"  : 'True',
      "MONOCLINIC" : 'alpha == gamma == 90',
      "ORTHORHOMBIC" : 'alpha == beta == gamma == 90',
      "TETRAGONAL" : 'a == b and alpha == beta == gamma == 90',
      "TRIGONAL"   : 'a == b == c and alpha == beta == gamma or ' +
                     'a == b and alpha == beta == 90 and gamma == 120',
      "HEXAGONAL"  : 'a == b and alpha == beta == 90 and gamma == 120',
      "CUBIC"      : 'a == b == c and alpha == beta == gamma == 90',
    }
    rule = crystal_system_rules[spacegroup.crystal_system]
    return eval(rule, locals())

# End of isSpaceGroupLatPar

# Helper class intended for this module only:
class Position2Tuple:
    """Create callable object that converts fractional coordinates to
    a tuple of integers with given precision.  For presision close to zero
    it will return a tuples of double.

    Data members:

    eps -- cutoff for equivalent coordinates.  When two coordiantes map to the
           same tuple, they are closer than eps.
    """
    def __init__(self, eps=epsilon):
        """Initialize Position2Tuple

        eps -- cutoff for equivalent coordinates
        """
        # ensure self.eps has exact machine representation
        self.eps = eps + 1.0
        self.eps = self.eps - 1.0
        # no conversions for very small eps
        if self.eps == 0.0 or 1.0/self.eps > sys.maxint:
            self.eps = 0.0
        return

    def __call__(self, xyz):
        """Convert array of fractional coordinates to a tuple.

        xyz -- fractional coordinates

        Return a tuple of 3 numbers.
        """
        # no conversion case
        if self.eps == 0.0:
            tpl = tuple(xyz % 1.0)
            return tpl
        # here we convert to integer
        tpl = tuple( [int((xi - numpy.floor(xi))/self.eps) for xi in xyz] )
        return tpl

# End of class Position2Tuple

def positionDifference(xyz0, xyz1):
    """Smallest difference between two coordinates in periodic lattice.

    xyz0, xyz1  -- fractional coordinates

    Return dxyz, a numpy.array dxyz with 0 <= dxyz <= 0.5.
    """
    dxyz = numpy.array(xyz0) - xyz1
    # map differences to [0,0.5]
    dxyz = dxyz - numpy.floor(dxyz)
    mask = (dxyz > 0.5)
    dxyz[mask] = 1.0 - dxyz[mask]
    return dxyz

# End of positionDifference

def nearestSiteIndex(sites, xyz):
    """Index of the nearest site to a specified position.

    sites -- list of coordinates or a 2-dimensional numpy.array
    xyz   -- single position

    Return integer.
    """
    # we use box distance to be consistent with Position2Tuple conversion
    dbox = positionDifference(sites, xyz).max(axis=1)
    nearindex = numpy.argmin(dbox)
    return nearindex

# End of nearestSiteIndex

def equalPositions(xyz0, xyz1, eps=epsilon):
    """Equality of two coordinates with optional tolerance.

    xyz0, xyz1 -- fractional coordinates
    eps        -- tolerance for equality of coordinates

    Return bool.
    """
    # we use box distance to be consistent with Position2Tuple conversion
    dxyz = positionDifference(xyz0, xyz1)
    return numpy.all(dxyz <= eps)

# End of equalPositions

def expandPosition(spacegroup, xyz, eps=epsilon):
    """Obtain unique equivalent positions and corresponding operations.

    spacegroup -- instance of SpaceGroup
    xyz        -- expanded position
    eps        -- cutoff for equal positions

    Return a tuple with (list of unique equivalent positions, nested
    list of SpaceGroups.SymOp instances, site multiplicity)
    """
    pos2tuple = Position2Tuple(eps)
    positions = []
    site_symops = {}    # position tuples with [related symops]
    for symop in spacegroup.iter_symops():
        pos = symop(xyz)
        tpl = pos2tuple(pos)
        if not tpl in site_symops:
            pos_is_new = True
            site_symops[tpl] = []
            # double check if there is any position nearby
            if positions:
                nearpos = positions[nearestSiteIndex(positions, pos)]
                # is it an equivalent position?
                if equalPositions(nearpos, pos, eps):
                    # tpl should map to the same list as nearpos
                    site_symops[tpl] = site_symops[ pos2tuple(nearpos) ]
                    pos_is_new = False
            if pos_is_new:  positions.append(pos)
        # here tpl is inside site_symops
        site_symops[tpl].append(symop)
    # pos_symops is nested list of symops associated with each position
    pos_symops = [ site_symops[pos2tuple(pos)] for pos in positions ]
    multiplicity = len(positions)
    return positions, pos_symops, multiplicity

# End of expandPosition

def expandAsymmetricUnit(spacegroup, asymunit, eps=0.0):
    """Expand positions in the asymmetric unit.

    spacegroup   -- instance of SpaceGroup
    asymunit     -- list of positions in asymmetric unit, it may
                    contain duplicates
    eps          -- cutoff for duplicate positions

    Return a tuple of (expandedunit, multiplicities).
    expandedunit -- list of equivalent positions per each site in asymunit
    multiplicities -- multiplicity of sites in asymunit
    """
    # By design Atom instances are not accepted as arguments to keep the
    # number of required imports low.
    expandedunit = []
    multiplicities = []
    for xyz in asymunit:
        eqsites, ignore, m = expandPosition(spacegroup, xyz, eps)
        expandedunit.append(eqsites)
        multiplicities.append(m)
    return expandedunit, multiplicities

# End of expandAsymmetricUnit

def nullSpace(A):
    """Null space of matrix A.
    """
    from numpy import linalg
    u, s, v = linalg.svd(A)
    # s may have smaller dimension than v
    mask = numpy.array([True]*numpy.shape(v)[0])
    mask[s > epsilon] = False
    null_space = numpy.compress(mask, v, axis=0)
    return null_space

# End of nullSpace

class GeneratorSite:
    """Storage of data related to a generator positions

    Data members:
        xyz          -- fractional coordinates of generator position
        eps          -- cutoff for equal positions
        eqxyz        -- list of equivalent positions
        symops       -- nested list of operations per each eqxyz
        multiplicity -- generator site multiplicity
        invariants   -- list of invariant operations for generator site
        null_space   -- null space for all possible differences of rotational
                        matrices from invariant operations.
        variables    -- list of (xyz symbol, value) pairs
    """

    def __init__(self, spacegroup, xyz, eps=epsilon):
        """Initialize GeneratorSite.

        spacegroup -- instance of SpaceGroup
        xyz        -- generating site.  When xyz is close to special
                      position self.xyz will be adjusted.
        eps        -- cutoff for equal positions
        """
        # just declare the variables
        self.xyz = None
        self.eps = eps
        self.eqxyz = None
        self.symops = None
        self.multiplicity = None
        self.invariants = []
        self.null_space = None
        self.variables = []
        # fill in the values
        self.xyz = xyz
        sites, ops, mult = expandPosition(spacegroup, xyz, eps)
        # shift self.xyz exactly to the special position
        if mult > 1:
            xyzdups = numpy.array([op(xyz) for op in ops[0]])
            dxyz = xyzdups - xyz
            dxyz = numpy.mean(dxyz - dxyz.round(), axis=0)
            # recalculate if needed
            if numpy.any(dxyz != 0.0):
                self.xyz = xyz + dxyz
                sites, ops, mult = expandPosition(spacegroup, self.xyz, eps)
        # self.xyz, sites, ops are all adjusted here
        self.eqxyz = sites
        self.symops = ops
        self.multiplicity = mult
        # invariant operations are always first in self.symop
        self.invariants = self.symops[0]
        self._findNullSpace()
        self._findVariables()
        return

    def _findNullSpace(self):
        """Calculate self.null_space from self.invariants.
        Try to represent self.null_space using small integers.
        """
        R0 = self.invariants[0].R
        Rdiff = [ (symop.R - R0) for symop in self.invariants ]
        Rdiff = numpy.concatenate(Rdiff, 0)
        self.null_space = nullSpace(Rdiff)
        if self.null_space.size == 0:   return
        # reverse sort rows of null_space rows by absolute value
        key = tuple(numpy.fabs(numpy.transpose(self.null_space))[::-1])
        order = numpy.lexsort(key)
        self.null_space = self.null_space[order[::-1]]
        # rationalize by the smallest element larger than cutoff
        cutoff = 1.0/32
        for i in range(len(self.null_space)):
            row = self.null_space[i]
            small = numpy.fabs(row[numpy.fabs(row) > cutoff]).min()
            signedsmall = row[numpy.fabs(row) == small][0]
            self.null_space[i] = self.null_space[i] / signedsmall
        return

    def _findVariables(self):
        """Find necessary variables and their values for expressing self.xyz
        """
        usedsymbol = {}
        # variable values depend on offset of self.xyz
        txyz = self.xyz
        # define txyz such that most of its elements are zero
        for nvec in self.null_space:
            idx = numpy.where(numpy.fabs(nvec) >= epsilon)[0][0]
            varvalue = txyz[idx]/nvec[idx]
            txyz = txyz - varvalue*nvec
            # determine variable name
            vname = [s for s in "xyz"[idx:] if not s in usedsymbol][0]
            self.variables.append( (vname, varvalue) )
            usedsymbol[vname] = True
        return

    def positionFormula(self, pos, xyzsymbols=("x","y","z")):
        """Formula of equivalent position with respect to generator site

        pos        -- fractional coordinates of possibly equivalent site
        xyzsymbols -- symbols for parametrized coordinates

        Return tuple of (xformula, yformula, zformula) formulas or empty tuple
        when pos is not equivalent to generator.  Formulas are formatted as
        "[[-][%g*]{x|y|z}] [{+|-}%g]", for example "z +0.5", "0.25", where
        <space> is required between parameter and constant part.
        """
        # find pos in eqxyz
        idx = nearestSiteIndex(self.eqxyz, pos)
        eqpos = self.eqxyz[idx]
        if not equalPositions(eqpos, pos, self.eps):    return ( )
        # any rotation matrix should do fine
        R = self.symops[idx][0].R
        nsrotated = numpy.dot(self.null_space, numpy.transpose(R))
        # build formulas using eqpos
        # find offset
        teqpos = numpy.array(eqpos)
        for nvec, (vname, varvalue) in zip(nsrotated, self.variables):
            teqpos -= nvec * varvalue
        # map varnames to xyzsymbols
        name2sym = dict( zip(("x", "y", "z"), xyzsymbols) )
        xyzformula = 3*[""]
        for nvec, (vname, varvalue) in zip(nsrotated, self.variables):
            for i in range(3):
                if abs(nvec[i]) < epsilon:  continue
                xyzformula[i] += "%+g*%s " % (nvec[i], name2sym[vname])
        # add constant offset teqpos to all formulas
        for i in range(3):
            if xyzformula[i] and teqpos[i] < epsilon: continue
            xyzformula[i] += "%+g" % teqpos[i]
        # reduce unnecessary +1* and -1*
        xyzformula = [ re.sub('^[+]1[*]|(?<=[+-])1[*]', '', f).strip()
                       for f in xyzformula ]
        return tuple(xyzformula)

# End of class GeneratorSite

def positionConstraints(spacegroup, positions, eps=epsilon, xyzsymbols=None):
    """Obtain symmetry constraints for specified positions.

    Procedure starts with initial list of parameter symbols of
    [ ["x0","y0","z0"], ["x1","y1","z1"], ... ] which is reduced
    to contain only independent parameters.

    spacegroup -- instance of SpaceGroup
    positions  -- list of all positions to be constrained
    eps        -- cutoff for equivalent positions
    xyzsymbols -- custom symbols for "x", "y", "z" per every coordinate

    Return a tuple of (poseqns, variables), where

    poseqns    -- list of coordinate formulas.  Formulas are formatted
                  as [[-]{x|y|z}%i] [{+|-}%g], for example: "x0", "-x3",
                  "z7 +0.5", "0.25".  Space before constant is required.
    variables  -- list of (xyz symbol, value) pairs
    """
    # positions returned by expandAsymmetricUnit have 3 dimensions
    positions = numpy.array(positions)
    if positions.ndim == 3:     positions = numpy.concatenate(positions, 0)
    numpos = len(positions)
    # check xyzsymbols
    if xyzsymbols is None:
        xyzsymbols = [ smbl+str(i) for i in range(numpos) for smbl in "xyz" ]
    if len(xyzsymbols) < positions.size:
        emsg = "Not enough xyz symbols for %i coordinates" % positions.size
        raise RuntimeError, emsg
    # we should be fine here
    poseqns = numpos*[[]]
    variables = []
    independent = dict.fromkeys(range(numpos))
    for genidx in range(numpos):
        if not genidx in independent:   continue
        # it is a generator
        genpos = positions[genidx]
        gen = GeneratorSite(spacegroup, genpos, eps)
        # append new variable if there are any
        gensymbols = xyzsymbols[3*genidx : 3*(genidx+1)]
        for k, v in gen.variables:
            smbl = gensymbols["xyz".index(k)]
            variables.append( (smbl, v) )
        # search for equivalents inside indies
        indies = independent.keys()
        indies.sort()
        for indidx in indies:
            indpos = positions[indidx]
            formula = gen.positionFormula(indpos, gensymbols)
            # formula is empty when indidx is independent
            if not formula:  continue
            # indidx is dependent here
            del independent[indidx]
            poseqns[indidx] = formula
    # all done here
    return poseqns, variables

# basic demonstration
if __name__ == "__main__":
    from SpaceGroups import sg100
    site = [.125, .625, .13]
    g = GeneratorSite(sg100, site)
    fm100 = g.positionFormula(site)
    print "g = GeneratorSite(sg100, %r)" % site
    print "g.positionFormula(%r) = %s" % (site, fm100)
    print "g.variables =", g.variables

# End of file

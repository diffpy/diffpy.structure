"""Symmetry utility functions such as expansion of asymmetric unit
"""

# version
__id__ = '$Id$'

import sys
import re
import numpy

def isSpaceGroupLatPar(spacegroup, a, b, c, alpha, beta, gamma):
    """Check if space group allows passed lattice parameters

    spacegroup -- instance of SpaceGroup
    a, b, c, alpha, beta, gamma -- lattice parameters

    return bool
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

# helper class to be used only inside this module
class Position2Tuple:
    """Create callable object that converts fractional coordinates to
    a tuple of integers with given precision.  For very high presision it
    returns tuples of doubles.

    Data members:
    eps      -- cutoff for duplicate coordinates.  Positions closer than eps
                yield intersecting pairs of tuples.
    """
    def __init__(self, eps):
        """Initialize Position2Tuple given the 

        eps  -- cutoff for duplicate coordinates
        """
        # ensure self.eps has exact machine representation
        self.eps = eps + 1.0
        self.eps = self.eps - 1.0
        # no conversions for very small eps
        if self.eps == 0.0 or 1.0/self.eps > sys.maxint:
            self.eps = 0.0
        return

    def __call__(self, xyz):
        """Low and high tuple hashes of fractional coordinates. 
        Coordinates xyz are mapped to 0.0 <= xyz < 1.0.

        xyz -- fractional coordinates

        Returns a tuple of integer equivalents.
        """
        # no conversion case
        if self.eps == 0.0:
            tpl = tuple(xyz % 1.0)
            return tpl
        # here we convert to integer
        tpl = tuple( [int((xi - numpy.floor(xi))/self.eps) for xi in xyz] )
        return tpl

# End of Position2Tuple

def positionDifference(xyz0, xyz1):
    """Smallest difference between 2 sites in periodic lattice.

    xyz0, xyz1  -- fractional coordinates

    return dxyz -- numpy.array where 0 <= dxyz <= 0.5
    """
    dxyz = numpy.array(xyz0) - xyz1
    # map differences to [0,0.5]
    dxyz = dxyz - numpy.floor(dxyz)
    mask = (dxyz > 0.5)
    dxyz[mask] = 1.0 - dxyz[mask]
    return dxyz

# End of positionDifference

def nearestSiteIndex(sites, xyz):
    """Index of the nearest site to a specified position

    sites -- list of positions or 2-dimensional numpy.array
    xyz   -- single position

    return integer
    """
    dsites = positionDifference(sites, len(sites)*[xyz])
    nearindex = numpy.argmin(numpy.sum(dsites,1))
    return nearindex

# End of nearestSiteIndex

def equalPositions(xyz0, xyz1, eps=0.0):
    """Equivalence of two positions with optional tolerance.

    xyz0, xyz1 -- fractional coordinates
    eps        -- tolerance on coordinate difference

    return bool
    """
    dxyz = positionDifference(xyz0, xyz1)
    return numpy.all(dxyz <= eps)

# End of equalPositions

def expandPosition(spacegroup, xyz, eps=0.0):
    """Obtain unique equivalent positions and corresponding operations.

    spacegroup -- instance of SpaceGroup
    xyz        -- original position
    eps        -- cutoff for duplicate positions

    returns a tuple with (list of unique equivalent positions, nested
    list of SpaceGroups.SymOp instances, site multiplicity)
    """
    pos2tuple = Position2Tuple(eps)
    positions = []
    site_symops = {}    # hashed position : related symops
    for symop in spacegroup.iter_symops():
        pos = symop(xyz)
        tpl = pos2tuple(pos)
        if not tpl in site_symops:
            pos_is_new = True
            site_symops[tpl] = []
            # but double check if there is any position nearby
            if positions:
                nearindex = nearestSiteIndex(positions, pos)
                # is it an equivalent position?
                if equalPositions(positions[nearindex], pos):
                    # tpl should map to the same list as neartpl
                    neartpl = pos2tuple(positions[nearindex])
                    site_symops[tpl] = site_symops[neartpl]
                    pos_is_new = False
            if pos_is_new:  positions.append(pos)
        # here tpl is inside site_symops
        site_symops[tpl].append(symop)
    # pos_symops contains symops associated with each position
    pos_symops = [ site_symops[pos2tuple(pos)] for pos in positions ]
    multiplicity = len(pos_symops[0])
    return positions, pos_symops, multiplicity

# End of expandPosition

def expandAsymmetricUnit(spacegroup, asymunit, eps=0.0):
    """Obtain unique equivalent positions and corresponding operations.

    spacegroup  -- instance of SpaceGroup
    asymunit    -- list of positions in asymmetric unit, it may
                   contain duplicates
    eps         -- cutoff for duplicate positions

    returns a nested list of equivalent positions, per each site in asymunit
    """
    # By design Atom instances are not accepted so that the number
    # of required modules is low.
    expanded = []
    for xyz in asymunit:
        eqsites = expandPosition(spacegroup, xyz, eps)[0]
        expanded.append(eqsites)
    return expanded

# End of expandAsymmetricUnit

# machine precision
epsilon = 1.0e-15

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

class GeneratorSite:
    """Storage of data related to a generator positions

    Data members:
        xyz          -- fractional coordinates of generator position
        eqxyz        -- list of equivalent positions
        symops       -- list of equivalent positions
        multiplicity -- generator site multiplicity (Wyckoff number)
        invariants   -- list of invariant operations for xyz
        null_space   -- numpy.array null space for differences of rotational
                        matrices from invariant operations.
        linked       -- list for indices of linked positions and related
                        transformations - to be filled by GeneratorSite user
    """

    def __init__(self, spacegroup, xyz, eps=0.0):
        """Initialize GeneratorSite

        spacegroup -- instance of SpaceGroup
        xyz        -- fractional coordinates of site nearest to [0,0,0]
        eps        -- cutoff for duplicate positions
        """
        self.xyz = None
        self.eps = eps
        self.eqxyz = None
        self.symops = None
        self.multiplicity = None
        self.invariants = []
        self.null_space = None
        self.variables = {}
        self.varnames = []
        self.linked = []
        self.xyz = numpy.array(xyz)
        # calculate multiplicity and invariants
        sites, ops, mult = expandPosition(spacegroup, xyz, eps)
        self.eqxyz = sites
        self.symops = ops
        self.multiplicity = mult
        # invariant operations come always first in self.symop
        self.invariants = self.symops[0]
        self._findNullSpace()
        self._findVariables()
        return

    def _findNullSpace(self):
        """Calculate null_space from self.invariants
        """
        R0 = self.invariants[0].R
        Rdiff = [ (symop.R - R0) for symop in self.invariants ]
        Rdiff = numpy.concatenate(Rdiff, 0)
        self.null_space = nullSpace(Rdiff)
        # reverse sort null_space rows by absolute value
        key = tuple(numpy.fabs(numpy.transpose(self.null_space))[::-1])
        order = numpy.lexsort(key)
        self.null_space = self.null_space[order[::-1]]
        # rationalize by the smallest, nonzero element
        cutoff = 1.0/32
        for i in range(len(self.null_space)):
            row = self.null_space[i]
            small = min( numpy.fabs(row[numpy.fabs(row) > cutoff]) )
            signedsmall = row[numpy.fabs(row) == small][0]
            self.null_space[i] = self.null_space[i] / signedsmall
        return

    def _findVariables(self):
        """Determine txyz from null_space
        """
        nsrank = numpy.rank(self.null_space)
        # xyz offset txyz cannot be expressed using null_space vectors
        self.txyz = self.xyz
        for nvec in self.null_space:
            idx = numpy.where(numpy.fabs(nvec) >= epsilon)[0][0]
            projection = self.txyz[idx]/nvec[idx]
            self.txyz = self.txyz - projection*nvec
            # determine variable name
            vname = "xyz"[idx]
            self.variables[vname] = projection
        self.varnames = self.variables.keys()
        self.varnames.sort()
        return

    def positionFormula(self, pos, xyzsymbols=("x","y","z")):
        """Formula of equivalent position with respect to generator position

        pos        -- fractional coordinates of possibly equivalent site
        xyzsymbols -- symbols for generator position

        return tuple of (xformula, yformula, zformula) formulas or empty tuple when
        pos is not equivalent to generator.
        """
        # find pos in eqxyz
        idx = nearestSiteIndex(self.eqxyz, pos)
        if not equalPositions(self.eqxyz[idx], pos):    return ( )
        # any rotation matrix should do fine
        R = self.symops[idx][0].R
        nsrot = numpy.dot(self.null_space, numpy.transpose(R))
        tpos = numpy.array(pos)
        for nvec, vname in zip(list(nsrot), self.varnames):
            tpos -= nvec*self.variables[vname] 
        # build formulas
        # map varnames to xyzsymbols
        name2sym = dict( zip(("x", "y", "z"), xyzsymbols) )
        xyzformula = 3*[""]
        for nvec, vname in zip(list(nsrot), self.varnames):
            for i in range(3):
                if abs(nvec[i]) < epsilon:  continue
                xyzformula[i] += "%+g*%s " % (nvec[i], name2sym[vname])
        # add constant offset tpos to the formulas
        for i in range(3):
            if xyzformula[i] and tpos[i] < epsilon: continue
            xyzformula[i] += "%+g" % tpos[i]
        # reduce +1* and -1*
        xyzformula = [ re.sub('^[+]1[*]|(?<=-)1[*]', '', f).strip()
                        for f in xyzformula ]
        return tuple(xyzformula)

# End of GeneratorSite

def positionConstraints(spacegroup, positions, eps=0.0):
    """Obtain symmetry constraints for specified positions.

    Procedure starts with initial list of constraints of
    [ ["x0","y0","zo"], ["x1","y1","z1"], ... ] and it reduces
    this list to contain only independent variables.

    spacegroup -- instance of SpaceGroup
    positions  -- list of all positions in to be constrained
    eps        -- cutoff for equivalent positions

    returns a tuple of (poseqns, variables) where

    poseqns    -- nested list of coordinate formulas.  Formulas are formatted
                  as [{+|-}{x|y|z}i] [{+|-}c], for example: "x0", "-x3",
                  "z7 +0.5", "0.25".
    variables  -- dictionary of variable names and values
    """
# FIXME
    # generators maps index of independent position to GeneratorSite
    generators = {}
    numpos = len(positions)
    independent = dict.fromkeys(range(numpos))
    for genidx in range(numpos):
        if not genidx in independent:   continue
        # it is a generator
        del independent[genidx]
        generators[genidx] = gen = GeneratorSite()
        gen.xyz = positions[genidx]
        eqpos, symops, mults = expandPosition(spacegroup, gen.xyz, eps)
        # find where in eqpos is gen.xyz
        geneqposidx = eqpos.index(gen.xyz)
        gen.multiplicity = mults[geneqposidx]
        gen.invariants = symops[geneqposidx]
        # search for equivalents inside indies
        indies = independent.keys()
        indies.sort()
        for indidx in indies:
            indpos = positions[indidx]
            for eqidx in range(len(eqpos)):
                if not equalPositions(indpos, eqpos[eqidx]):    continue
                # indpos is linked to generator via any of symops[eqidx]
                gen.linked.append( (indidx, symops[eqidx][0]) )
                del independent[indidx]
                # no need to search eqpos anymore
                break
    # all generators should be found here


# basic demonstration
if __name__ == "__main__":
    from SpaceGroups import sg100
    site = [.125, .625, .13]
    g = GeneratorSite(sg100, site)
    print "g = GeneratorSite(sg100, %r)" % site
    print "g.positionFormula(%r) = %s" % (site, g.positionFormula(site))
    print "g.variables =", g.positionFormula(site)

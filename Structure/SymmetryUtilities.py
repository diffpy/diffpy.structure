"""Symmetry utility functions such as expansion of asymmetric unit
"""

import sys
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

def expandPosition(spacegroup, xyz, eps=1E-8):
    """Obtain unique equivalent positions and corresponding operations.

    spacegroup -- instance of SpaceGroup
    xyz        -- original position
    eps        -- cutoff for duplicate positions

    returns a tuple with (list of unique equivalent positions, nested
    list of SpaceGroups.SymOp instances, list of site multiplicities)
    """
    pos2tuple = Position2Tuple(eps)
    positions = []
    site_symops = { }   # hashed position : related symops
    for symop in spacegroup.iter_symops():
        pos = symop(xyz)
        tpl = pos2tuple(pos)
        if not tpl in site_symops:
            pos_is_new = True
            site_symops[tpl] = []
            # but double check if there is any position nearby
            if positions:
                dpos = numpy.array([p-pos for p in positions])
                dpos = dpos - numpy.floor(dpos)
                dpos[dpos > 0.5] = 1.0 - dpos[dpos > 0.5]
                nearindex = numpy.argmin(sum(dpos,1))
                dnear = dpos[nearindex]
                # there is an equivalent position
                if numpy.all(dnear < eps):
                    # tpl should map to existing list
                    neartpl = pos2tuple(positions[nearindex])
                    site_symops[tpl] = site_symops[neartpl]
                    pos_is_new = False
            if pos_is_new:  positions.append(pos)
        # here tpl is inside site_symops
        site_symops[tpl].append(symop)
    # pos_symops contains symops associated with each position
    pos_symops = [ site_symops[pos2tuple(pos)] for pos in positions ]
    multiplicities = [ len(lst) for lst in pos_symops ]
    return positions, pos_symops, multiplicities

# End of expandPosition

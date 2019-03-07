#!/usr/bin/env python
##############################################################################
#
# diffpy.structure  by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2008 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Chris Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################

"""Utilities for making shapes."""

def findCenter(S):
    """Find the approximate center atom of a structure.

    The center of the structure is the atom closest to (0.5, 0.5, 0.5)

    Returns the index of the atom.
    """
    best = -1
    bestd = len(S)
    center = [0.5, 0.5, 0.5] # the cannonical center

    for i in range(len(S)):
        d = S.lattice.dist(S[i].xyz, center)
        if d < bestd:
            bestd = d
            best = i

    return best

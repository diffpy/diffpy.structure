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
    """Find the approximate center `Atom` of a `Structure`.

    The center of the `Structure` is the `Atom` closest to ``(0.5, 0.5, 0.5)``.

    Parameters
    ----------
    S : Structure
        A `Structure` instance.

    Returns
    -------
    int
        The index of the center `Atom`.
    """
    best = -1
    bestd = len(S)
    center = [0.5, 0.5, 0.5]  # the cannonical center

    for i in range(len(S)):
        d = S.lattice.dist(S[i].xyz, center)
        if d < bestd:
            bestd = d
            best = i

    return best

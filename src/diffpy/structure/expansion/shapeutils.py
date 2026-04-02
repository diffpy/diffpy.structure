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
from diffpy.utils._deprecator import build_deprecation_message, deprecated

base = "diffpy.structure"
removal_version = "4.0.0"
findCenter_deprecation_msg = build_deprecation_message(
    base,
    "findCenter",
    "find_center",
    removal_version,
)
"""Utilities for making shapes."""


@deprecated(findCenter_deprecation_msg)
def findCenter(S):
    """This function has been deprecated and will be removed in version
    4.0.0.

    Please use diffpy.structure.find_center instead.
    """
    return find_center(S)


def find_center(S):
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
    center = [0.5, 0.5, 0.5]  # the canonical center

    for i in range(len(S)):
        d = S.lattice.dist(S[i].xyz, center)
        if d < bestd:
            bestd = d
            best = i

    return best

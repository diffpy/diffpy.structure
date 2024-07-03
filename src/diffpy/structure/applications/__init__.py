#!/usr/bin/env python
##############################################################################
#
# diffpy.structure  Complex Modeling Initiative
#                   (c) 2019 Brookhaven Science Associates,
#                   Brookhaven National Laboratory.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################

"""
Import support for renamed subpackage diffpy.structure.applications.

This module is deprecated and will be removed in version 3.1.
"""

# TODO remove this module in version 3.1.

from warnings import warn

warn(
    "Module 'diffpy.structure.applications' is deprecated.  " "Import 'diffpy.structure.apps' instead.",
    DeprecationWarning,
    stacklevel=2,
)

#!/usr/bin/env python
##############################################################################
#
# diffpy.structure  Complex Modeling Initiative
#                   (c) 2017 Brookhaven Science Associates,
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
Support import of old camel-case module names with DeprecationWarning.

The imported camel-case modules are aliases for the current module
instances.  Their `__name__` attributes are thus all in lower-case.

This module is deprecated and will be removed in the future.
"""


import sys

# install legacy import hooks
import diffpy.structure._legacy_importer

# replace this module with the new one
sys.modules['diffpy.Structure'] = diffpy.structure

# End of file

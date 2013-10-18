#!/usr/bin/env python
##############################################################################
#
# diffpy.Structure  by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2008 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################


"""Definition of __version__, __date__, __gitsha__ for diffpy.Structure.
"""


# obtain version information
from pkg_resources import get_distribution
__version__ = get_distribution('diffpy.Structure').version

# we assume that tag_date was used and __version__ ends in YYYYMMDD
__date__ = __version__[-8:-4] + '-' + \
           __version__[-4:-2] + '-' + __version__[-2:]

# GIT SHA hash is the second last component in the version string
__gitsha__ = __version__.rsplit('-', 2)[-2][1:]

# End of file

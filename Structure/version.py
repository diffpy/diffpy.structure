########################################################################
#
# Structure         by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2006 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
########################################################################

"""Definition of __version__ for Structure container package.
"""

__id__ = "$Id$"

try:
    from diffpy.version import __svnrevision__
except ImportError:
    __svnrevision__ = "?"

__version__ = "0.1." + __svnrevision__

# End of file

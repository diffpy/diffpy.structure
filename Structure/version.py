########################################################################
#
# <PackageName>     by DANSE Diffraction group
#                   Simon J.L. Billinge
#                   Michigan State University
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See COPYRIGHT.txt for copying and usage conditions.
# See LICENSE.txt for license information.
#
########################################################################

"""Definition of __version__ for Structure container package.
"""

# to update __version__ add hashmark below just before committing
##

__id__ = "$Id$"

__version__ = "0.1."

svnrevision = __id__.split()[2:3]
if svnrevision:     __version__ += svnrevision[0]
else:               __version__ += "?"

# End of file

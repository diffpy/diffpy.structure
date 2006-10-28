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

"""Definition of version for Structure container package.
"""

# to update version add hashmark below just before committing
##

__id__ = "$Id$"

version = "0.1."
svnrevision = __id__.split()[2:3]
if svnrevision:     version += svnrevision[0]
else:               version += "?"

# End of file

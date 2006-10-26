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

"""definition of PDFFitStructure class derived from Structure
"""

__id__ = "$Id$"

from Structure import Structure

##############################################################################
class PDFFitStructure(Structure):
    """PDFFitStructure --> Structure with extra pdffit member

    Data members:
        pdffit -- dictionary for storing following extra parameters from
                  PDFFit structure files:
                      'scale', 'delta', 'gamma', 'srat',
                      'rcut', 'spcgr', 'dcell', 'ncell'
    """

    def __init__(self, *args, **kwargs):
        Structure.__init__(self, *args, **kwargs)
        self.pdffit = {
            'scale'  :  1.0,
            'delta'  :  0.0,
            'gamma'  :  0.0,
            'srat'   :  1.0,
            'rcut'   :  0.0,
            'spcgr'  :  'P1',
            'dcell'  :  6*[0.0],
            'ncell'  :  [1, 1, 1, 0],
        }
        return

# End of PDFFitStructure

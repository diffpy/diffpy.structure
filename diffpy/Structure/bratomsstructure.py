#!/usr/bin/env python
##############################################################################
#
# diffpy.Structure  by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2006 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    C. L. Farrow
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################

"""definition of BRAtomsStructure class derived from Structure
"""

from diffpy.Structure.structure import Structure

##############################################################################
class BRAtomsStructure(Structure):
    """BRAtomsStructure --> Structure with extra information for use with Bruce
    Ravel's atoms program.

    Data members:
        bratoms -- dictionary for storing following extra parameters from
                   atoms .inp files.
                   'space', 'output', 'rmax', 'core', 'edge', 'shift',
                   'nitrogen', 'argon', 'krypton'
        see the following web site for descriptions.
   http://leonardo.phys.washington.edu/~ravel/software/doc/Atoms/Atoms/node7.html
    """

    def __init__(self, *args, **kwargs):
        Structure.__init__(self, *args, **kwargs)
        self.bratoms = {
            'space'     : "p1",
            'rmax'      :  6.0,
            'core'      : None,
            'edge'      : 'K',
            'shift'     : None,
            'output'    : None,
            'nitrogen'  : None,
            'argon'     : None,
            'krypton'   : None,
        }
        return

# End of BRAtomsStructure

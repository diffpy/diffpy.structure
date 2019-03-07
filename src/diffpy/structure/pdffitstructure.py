#!/usr/bin/env python
##############################################################################
#
# diffpy.Structure  by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2006 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################

"""definition of PDFFitStructure class derived from Structure
"""


from diffpy.Structure.structure import Structure

##############################################################################
class PDFFitStructure(Structure):
    """PDFFitStructure --> Structure with extra pdffit member

    Data members:
        pdffit -- dictionary for storing following extra parameters from
                  PDFFit structure files:
                      'scale', 'delta1', 'delta2', 'sratio',
                      'rcut', 'spcgr', 'dcell', 'ncell'
    """

    def __init__(self, *args, **kwargs):
        self.pdffit = {
            'scale' : 1.0,
            'delta1' : 0.0,
            'delta2' : 0.0,
            'sratio' : 1.0,
            'rcut' : 0.0,
            'spcgr' : 'P1',
            'spdiameter' : 0.0,
            'stepcut' : 0.0,
            'dcell' : 6*[0.0],
            'ncell' : [1, 1, 1, 0],
        }
        Structure.__init__(self, *args, **kwargs)
        return

    def read(self, filename, format='auto'):
        """Same as Structure.read, but update spcgr value in
        self.pdffit when parser can get spacegroup.

        Return instance of StructureParser used to load the data.
        See Structure.read() for more info.
        """
        p = Structure.read(self, filename, format)
        sg = getattr(p, 'spacegroup', None)
        if sg:  self.pdffit['spcgr'] = sg.short_name
        return p

    def readStr(self, s, format='auto'):
        """Same as Structure.readStr, but update spcgr value in
        self.pdffit when parser can get spacegroup.

        Return instance of StructureParser used to load the data.
        See Structure.readStr() for more info.
        """
        p = Structure.readStr(self, s, format)
        sg = getattr(p, 'spacegroup', None)
        if sg:  self.pdffit['spcgr'] = sg.short_name
        return p

# End of PDFFitStructure

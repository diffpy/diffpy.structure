#!/usr/bin/env python
##############################################################################
#
# diffpy.structure  by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2007 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################

"""Index of recognized structure formats, their IO capabilities and
associated modules where they are defined.  Plugins for new structure
formats need to be added to the parser_index dictionary in this module.
"""

parser_index = {

    # automatic format detection - tries all parsers one by one
    'auto' : {
        'module' : 'p_auto',
        'file_extension' : '',
        'file_pattern' : '*.*',
        'has_input' : True,
        'has_output' : False,
        },

    # CIF format
    'cif' : {
        'module' : 'p_cif',
        'file_extension' : '.cif',
        'file_pattern' : '*.cif',
        'has_input' : True,
        'has_output' : True,
        },

    # PDB format
    'pdb' : {
        'module' : 'p_pdb',
        'file_extension' : '.pdb',
        'file_pattern' : '*.pdb',
        'has_input' : True,
        'has_output' : True,
        },

    # Discus structure format
    'discus' : {
        'module' : 'p_discus',
        'file_extension' : '.stru',
        'file_pattern' : '*.stru|*.rstr',
        'has_input' : True,
        'has_output' : True,
        },

    # PDFfit structure format
    'pdffit' : {
        'module' : 'p_pdffit',
        'file_extension' : '.stru',
        'file_pattern' : '*.stru|*.rstr',
        'has_input' : True,
        'has_output' : True,
        },

    # standard xyz file
    'xyz' : {
        'module' : 'p_xyz',
        'file_extension' : '.xyz',
        'file_pattern' : '*.xyz',
        'has_input' : True,
        'has_output' : True,
        },

    # raw xyz file (element labels optional)
    'rawxyz' : {
        'module' : 'p_rawxyz',
        'file_extension' : '.xyz',
        'file_pattern' : '*.xyz',
        'has_input' : True,
        'has_output' : True,
        },

    # AtomEye extended configuration format
    'xcfg' : {
        'module' : 'p_xcfg',
        'file_extension' : '',
        'file_pattern' : '*.xcfg|*.eye|*.cfg',
        'has_input' : True,
        'has_output' : True,
        },
}

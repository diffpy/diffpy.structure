# -*- Makefile -*-
########################################################################
#
# Structure         by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2007 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See COPYRIGHT.txt for copying and usage conditions.
# See LICENSE.txt for license information.
#
########################################################################

PROJECT = diffpy
PACKAGE = Structure/Parsers

#--------------------------------------------------------------------------
#

all: export

#--------------------------------------------------------------------------
#
# export

EXPORT_PYTHON_MODULES = \
    P_auto.py           \
    P_cif.py            \
    P_pdb.py            \
    P_pdffit.py         \
    P_rawxyz.py         \
    P_xcfg.py           \
    P_xyz.py            \
    StructureParser.py  \
    import_helper.py	\
    parser_index.py	\
    __init__.py

export:: export-package-python-modules

# version
# $Id$

# End of file

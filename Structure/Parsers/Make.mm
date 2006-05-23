# -*- Makefile -*-
#
# <LicenseText>

PROJECT = Structure/Parsers
PACKAGE = Parsers

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
    __init__.py

export:: export-python-modules

# version
# $Id$

# End of file

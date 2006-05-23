# -*- Makefile -*-
#
# <LicenseText>

PROJECT = Structure
PACKAGE = Structure

# directory structure

BUILD_DIRS = \
    Parsers

OTHER_DIRS =

RECURSE_DIRS = $(BUILD_DIRS) $(OTHER_DIRS)

#--------------------------------------------------------------------------
#

all: export
	BLD_ACTION="all" $(MM) recurse

RECURSE_DIRS = Parsers

#--------------------------------------------------------------------------
#
# export

EXPORT_PYTHON_MODULES = \
    __init__.py         \
    PDFFitStructure.py  \
    PeriodicTable.py    \
    atom.py 		\
    lattice.py 		\
    structure.py 	\
    Parsers

export:: export-python-modules

unittests:
	PYTHONPATH=$(dir $(PWD) ):$$PYTHONPATH ./alltests.py

# version
# $Id$

# End of file

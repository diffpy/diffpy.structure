# -*- Makefile -*-
########################################################################
#
# Structure         by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   Michigan State University
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See COPYRIGHT.txt for copying and usage conditions.
# See LICENSE.txt for license information.
#
########################################################################

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

EXPORT_PYTHON_MODULES =     \
    __init__.py		    \
    PDFFitStructure.py	    \
    PeriodicTable.py	    \
    Atom.py		    \
    Lattice.py		    \
    Structure.py	    \
    SpaceGroups.py	    \
    SymmetryUtilities.py    \
    utils.py		    \
    StructureErrors.py	    \
    version.py		    \
    Parsers

export:: export-python-modules

unittests:
	PYTHONPATH=$(dir $(PWD) ):$$PYTHONPATH ./alltests.py

# version
# $Id$

# End of file

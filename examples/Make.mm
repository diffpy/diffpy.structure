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
PACKAGE = examples

#--------------------------------------------------------------------------
#

all: clean

release: clean
	cvs release .

update: clean
	cvs update .


# version
# $Id$

# End of file

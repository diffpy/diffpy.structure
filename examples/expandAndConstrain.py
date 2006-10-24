#!/usr/bin/env python

import sys
import os.path
import re
import numpy

from Structure.SymmetryUtilities import *
from Structure.SpaceGroups import GetSpaceGroup

# Get and process coordinates
print "Enter xyz coordinates (can be fractions): "
reply = raw_input('> ')
# split reply on commas or whitespace
words = re.split('[,\s]+', reply.strip())
# prepend "1.0*" to avoid integer division
xyz = [ eval("1.0*" + e) for e in words[:3] ]

# Get space group for expansion
print "Enter space group (symbol or number): "
reply = raw_input('> ')
# can it be converted to number?
try:
    sgid = int(reply)
except ValueError:
    sgid = reply.strip()
sg_expansion = GetSpaceGroup(sgid)
# get optional offset
sg_expoffset = [0.0, 0.0, 0.0]
print "Space group origin offset (<Enter> if none): "
reply = raw_input('> ').strip()
if reply:
    # split reply on commas or whitespace
    words = re.split('[,\s]+', reply)
    # prepend "1.0*" to avoid integer division
    sg_expoffset = [ eval("1.0*" + e) for e in words[:3] ]

# Expand coordinate
xyz_expanded, ignore, mltp = expandPosition(sg_expansion, xyz, sg_expoffset)
print
print "Site multiplicity:", mltp
print "Expanded positions:"
for eqxyz in xyz_expanded:
    print "   ", str(eqxyz)
print

# Get space group for constrainment
print "Enter constraining space group: "
reply = raw_input('> ')
# can it be converted to number?
try:
    sgid = int(reply)
except ValueError:
    sgid = reply.strip()
sg_constrainment = GetSpaceGroup(sgid)
# get optional offset
sg_consoffset = [0.0, 0.0, 0.0]
print "Space group origin offset (<Enter> if none): "
reply = raw_input('> ').strip()
if reply:
    # split reply on commas or whitespace
    words = re.split('[,\s]+', reply)
    # prepend "1.0*" to avoid integer division
    sg_consoffset = [ eval("1.0*" + e) for e in words[:3] ]

# Generate constraint equations
print repr(sg_consoffset)
eqns, pars = positionConstraints(sg_constrainment, xyz_expanded, sg_consoffset)
print
print "Position formulas:"
for eqxyz in eqns:
    print "   ", str(eqxyz)
print "Parameters:"
for symbol, value in pars:
    print "   ", symbol, value

# PDFFIT2 macros
print
print "Display as PDFFIT macros? (y/n)"
reply = raw_input('> ')
if reply[:1].lower() != "y":
    sys.exit()

# check if string is float:
def isfloat(s):
    try:
        x = float(s)
        return True
    except ValueError:
        return False
# End of isfloat

# generate symbol for every possible coordinate
numcoordinates = 3*len(xyz_expanded)
firstpar = 11
parsymbols = [ "@"+str(i) for i in range(firstpar, firstpar+numcoordinates) ]

# create constraints using these symbols
eqns, pars = positionConstraints(sg_constrainment, xyz_expanded, sg_consoffset,
        xyzsymbols=parsymbols)

# print constrain commands
print 78*"#"
siteindex = 0
for (eqx, eqy, eqz) in eqns:
    siteindex += 1
    for varname, formula in [ ("x",eqx), ("y",eqy), ("z",eqz) ]:
        if isfloat(formula):    continue
        # here formula contains parameter
        print "constrain(%s(%i), %r)" % (varname, siteindex, formula)

# print setpar commands
for symbol, value in pars:
    print "setpar(%s, %s)" % ( symbol.lstrip('@'), value )

if not pars:
    print "# special fixed positions, no constraints"

print 78*"#"

# End of file

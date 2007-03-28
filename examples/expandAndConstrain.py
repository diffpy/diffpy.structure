#!/usr/bin/env python

import sys
import os.path
import re
import numpy

from diffpy.Structure.SymmetryUtilities import *
from diffpy.Structure.SpaceGroups import GetSpaceGroup

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
symcon = SymmetryConstraints(sg_constrainment, xyz_expanded, sgoffset=sg_consoffset)
print
print "Position formulas:"
for eq in symcon.positionFormulas():
    print "    (%r, %r, %r)" % (eq["x"], eq["y"], eq["z"])
print "Parameters:"
for symbol, value in symcon.pospars:
    print "   ", symbol, value

# PDFFIT2 macros
print
print "Display as PDFFIT macros? (y/n)"
reply = raw_input('> ')
if reply[:1].lower() != "y":
    sys.exit()

# generate symbol for every possible coordinate
numcoordinates = len(symcon.pospars)
firstpar = 11
parsymbols = [ "@"+str(i) for i in range(firstpar, firstpar+numcoordinates) ]

# create constraints using these symbols
eqns = symcon.positionFormulasPruned(parsymbols)

# print constrain commands
print 78*"#"
siteindex = 0
for eq in eqns:
    sortedkeys = eq.keys()
    sortedkeys.sort()
    siteindex += 1
    for smbl in sortedkeys:
        print "constrain(%s(%i), %r)" % (smbl, siteindex, eq[smbl])

# print setpar commands
for symbol, value in zip(parsymbols, symcon.posparValues()):
    print "setpar(%s, %s)" % ( symbol.lstrip('@'), value )

if not symcon.pospars:
    print "# special fixed positions, no constraints"

print 78*"#"

# End of file

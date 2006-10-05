"""This module includes functions to expand the cell structure to unitcell according to specified space group using cctbx package. Three functions are defined. These functions return a structure containing the  resulting structure and a string of constraints.
"""
import copy
import numpy
import Structure.structure 
from Structure.structure import Structure
from Structure.lattice import Lattice
from Structure.atom import Atom 

################################################################################

def expandStructureConstraints(S, clists):
    """Expand the structure S and its constraint equations.

    S       --  A structure to expand
    clists   --  List of constraint equation tuples. The list contains one tuple
                for each atom in the structure. The tuples are formatted as
                (x, y, z, u11, u22, u33, u12, u23, u13, occ). Each atom must
                have a tuple. Missing values in the tuple are ''.

    Returns a new (Structure, clists) tuple
    """
    if not isinstance(S, Structure):
        message = "Expected instance of Structure, got %s instead" % S
        raise SyntaxError, message
    if len(clists) != len(S):
        message =  "clist different length than structure (%i != %i)" \
                    % (len(clists), len(S))
        raise SyntaxError, message

    for clist in clists:
        if len(clist) != 10:
            message = "Constraint list has wrong length"
            raise SyntaxError, message



    Snew = expandStructure(S) 
    # Get the exact operations applied to get Snew
    symoplists = getEssentialOperations(S)

    # Create the new constraints list. 
    cnew = []
    # In order to update the space group equations, we have to apply the
    # rotation and translation operators of the essential symmetry operations to
    # the constraint equations. The u values and occ values don't change.
    anum = 0 # keep track of which of the new atoms we're dealing with
    for i in range(len(S)):
        symops = symoplists[i]
        # Get the constraint equations
        xyz = clists[i][:3]
        for j in range(len(symops)):
            cnew.append(['']*10)
            symop = symops[j]
            # Copy the u and occ constraint equations
            cnew[anum][3:] = clists[i][3:]
            # Apply the symmetry operations to the constraint equations
            t = symop.t
            R = symop.R
            newxyz = ['']*3

            # Loop over x, y, z in newxyz
            for k in range(3):
                # Loop over x, y, z in xyz (the rotation loop)
                for l in range(3):

                    if R[k][l] != 0 and xyz[l].strip():

                        # Format the addition to the pattern.
                        # Avoid patterns like '1.000*(@1)'
                        if R[k][l] == 1:
                            temp = xyz[l]
                        else:
                            temp = "%1.4f*(%s)" % (R[k][l], xyz[l])

                        # Add the addition if there is something to add it to
                        if newxyz[k]:
                            newxyz[k] = "%s+%s" % (newxyz[k], temp)
                        # Otherwise the current pattern is the new pattern
                        else:
                            newxyz[k] = temp

                # Apply the translation
                if t[k] != 0 and newxyz[k].strip():
                    if t[k] > 0:
                        newxyz[k] = "%s+%1.4f" % (newxyz[k], t[k])
                    else:
                        newxyz[k] = "%s%1.4f" % (newxyz[k], t[k])

            # Apply the new operation
            cnew[anum][:3] = newxyz[:3]
            anum += 1
                 

    return (Snew, cnew)

def getEssentialOperations(S):
    """Get the symmetry operations essential to unique lattice transforms.

    Returns a list of symmetry operation tuples. Each tuple corresponds to the
    essential operations for transforming a single atom.
    """
    zeros = numpy.zeros(3)
    spacegroup = S.lattice.spacegroup
    
    symoplists = []
    for atom in S:
        equivstr = []
        equivpos = []
        symops = []
        # Weed out duplicates and those xyzs with negative coordinates.
        # Is this correct?
        for symop in spacegroup.iter_symops():
            pos = symop(atom.xyz)
            if str(pos) not in equivstr and not (pos < zeros).any():
                equivpos.append(pos)
                equivstr.append(str(pos))
                symops.append(symop)
        symoplists.append(symops)

    return symoplists

def expandStructure(S):
    """Expand the structure based upon its space group.

    S       --  A structure to expand

    Returns a new Structure object.
    """
    symoplists = getEssentialOperations(S)
    # Want a new instance of whatever 'S' is. It may be a PDFStructure
    # object.
    Snew = copy.deepcopy(S)
    # Is there a faster way to clear the atoms?
    for i in range(len(S)): Snew.pop(-1)
    Snew.lattice.setSpaceGroup('P1')

    for i in range(len(S)):
        symops = symoplists[i]
        atom = S[i]
        #for symop in symops:
        for j in range(len(symops)):
            symop = symops[j]
            A = copy.deepcopy(atom)
            A.name = "%s_%i" % (atom.name, j)
            A.xyz = symop(atom.xyz)
            Snew.append(A)

    return Snew

################################################################################
    
def _printStructure(S):
    print S.title
    print "space group = %s" % S.lattice.spacegroup.short_name
    print len(S), "atoms"
    for atom in S:
        x, y, z = atom.xyz
        print "%s: %-2.4f %-2.4f %-2.4f" % (atom.name, x, y, z)
    return

def _testNi():
    L = Lattice(spacegroup="Fm-3m")
    S = Structure(lattice = L)
    S.title = "Ni"
    A0 = Atom(xyz = numpy.array([0,0,0]), element = "Ni", name = "Ni")
    S.append(A0)
    _printStructure(S)
    clist = ['']*10
    clist[0:3] = ["@1", "@2", "@3"]
    Sexp, clistsexp = expandStructureConstraints(S, [clist])
    _printStructure(Sexp)
    for l in clistsexp: print l
    return

def _testGaAs():
    L = Lattice(spacegroup="F-43m")
    S = Structure(lattice = L)
    S.title = "GaAs"
    clists = []
    A0 = Atom(xyz = numpy.array([0,0,0]), element = "Ga", name = "Ga")
    clist  = ['']*10
    clist[0:3] = ["@1", "@2", "@3"]
    clists.append(clist)
    A1 = Atom(xyz = numpy.array([0.5,0,0]), element = "As", name = "As")
    clist  = ['']*10
    clist[0:3] = ["@4", "@5", "@6"]
    clists.append(clist)
    S.append(A0)
    S.append(A1)
    _printStructure(S)
    Sexp = expandStructure(S)
    Sexp, clistsexp = expandStructureConstraints(S, clists)
    _printStructure(Sexp)
    for l in clistsexp: print l
    return

if __name__ == '__main__':
    _testNi()
    print "\n"
    #_testGaAs()

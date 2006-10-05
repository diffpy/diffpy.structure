"""This module includes functions to expand the cell structure to unitcell according to specified space group using cctbx package. Three functions are defined. These functions return a structure containing the  resulting structure and a string of constraints.
"""
import copy, types
import numpy
import Structure.structure 
from Structure.structure import Structure
from Structure.lattice import Lattice
from Structure.atom import Atom 

################################################################################

def expandStructureConstraints(S, clists):
    """Expand the structure S and its constraint equations.

    S       --  A structure to expand
    clists  --  List of constraint equation tuples. The list contains one tuple
                for each atom in the structure. The tuples are formatted as
                (x, y, z, u11, u22, u33, u12, u23, u13, occ). Each atom must
                have a tuple. Missing values in the tuple are ''.

    Returns a new (Structure, clists) tuple
    """
    if not isinstance(S, Structure):
        message = "Expected instance of Structure, got %s instead" % S
        raise TypeError, message
    if len(clists) != len(S):
        message =  "clist different length than structure (%i != %i)" \
                    % (len(clists), len(S))
        raise ValueError, message

    for clist in clists:
        if len(clist) != 10:
            message = "Constraint list has wrong length"
            raise ValueError, message

    # Get the exact operations applied to get Snew
    symoplists = getEssentialOperations(S)
    Snew = expandStructure(S, symoplists) 

    # Create the new constraints list. 
    cnew = []
    # In order to update the constraint equations, apply the rotation and
    # translation operators of the essential symmetry operations.  The u values
    # and occ values don't change.
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

def expandStructure(S, symoplists = None):
    """Expand the structure based upon its space group.

    S           --  A structure to expand
    symopilsts  --  A list of symmetry operations to be applied by each atom. If
                    symoplist is None (default), this is generated from the
                    space group.

    Returns a new Structure object.
    """
    if symoplists is None:
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

def getEssentialOperations(S):
    """Get the symmetry operations essential to unique lattice transforms.

    Returns a list of symmetry operation tuples. Each tuple corresponds to the
    essential operations for transforming a single atom.
    """
    zeros = numpy.zeros(3)
    ones = numpy.ones(3)
    spacegroup = S.lattice.spacegroup
    
    symoplists = []
    for atom in S:
        equivstr = []
        equivpos = []
        symops = []
        # Weed out duplicates and those xyzs that go outside the unit cell.
        for symop in spacegroup.iter_symops():
            pos = symop(atom.xyz)
            if str(pos) not in equivstr \
                and not (pos < zeros).any() \
                and not (pos >= ones).any():

                equivpos.append(pos)
                equivstr.append(str(pos))
                symops.append(symop)
        symoplists.append(symops)

    return symoplists

def makeSupercellConstraints(S, clists, l=1, m=1, n=1):
    """Make a supercell of structure S and its constraint equations.

    S       --  A structure to expand
    clists  --  List of constraint equation tuples. The list contains one tuple
                for each atom in the structure. The tuples are formatted as
                (x, y, z, u11, u22, u33, u12, u23, u13, occ). Each atom must
                have a tuple. Missing values in the tuple are ''.
    l       --  The positive integer multiplier along the a-axis (default 1)
    m       --  The positive integer multiplier along the b-axis (default 1)
    n       --  The positive integer multiplier along the c-axis (default 1)

    Returns a new (Structure, clists) tuple
    """
    if len(clists) != len(S):
        message =  "clist different length than structure (%i != %i)" \
                    % (len(clists), len(S))
        raise ValueError, message
    for clist in clists:
        if len(clist) != 10:
            message = "Constraint list has wrong length"
            raise ValueError, message

    Snew = makeSupercell(S, l, m, n)
    cnew = []
    lmn = [l, m, n]

    anum = 0 # the index of the new atom
    for a in range(len(S)):
        for i in range(l):
            for j in range(m):
                for k in range(n):
                    ijk = [i, j, k]
                    cnew.append(['']*10)
                    cnew[anum][3:] = clists[a][3:]
                    for h in range(3):
                        neweq = clists[a][h].strip()
                        if neweq: 
                            if ijk[h] != 0:
                                neweq = "%s+%i" % (neweq, ijk[h])
                            if lmn[h] != 1:
                                neweq = "(%s)/%i" % (neweq, lmn[h])
                        cnew[anum][h] = neweq
                    # for h
                    anum += 1
                # for k
            # for j
        # for i
    # for a

    return (Snew, cnew)

def makeSupercell(S, l=1, m=1, n=1):
    """Make a supercell out of a structure.

    S   --  The structure to expand
    l   --  The positive integer multiplier along the a-axis (default 1)
    m   --  The positive integer multiplier along the b-axis (default 1)
    n   --  The positive integer multiplier along the c-axis (default 1)

    The structure S should have a 'P1' supercell, but this is not enforced.

    Returns the new structure
    """
    if not isinstance(S, Structure):
        message = "Expected instance of Structure, got %s instead" % S
        raise TypeError, message
    message = "l, m, n must be positive integers"
    for i in (l, m, n):
        if i < 1:
            raise ValueError, message
        if not type(i) is types.IntType:
            raise ValueError, message

    L = S.lattice
    Snew = copy.deepcopy(S)
    Lnew = Lattice(l*L.a, m*L.b, n*L.c, L.alpha, L.beta, L.gamma, spacegroup = 'P1')
    # Is there a faster way to clear the atoms?
    for i in range(len(S)): Snew.pop(-1)
    Snew.lattice = Lnew

    for atom in S:
        dup = 0
        for i in range(l):
            for j in range(m):
                for k in range(n):
                    A = copy.deepcopy(atom)
                    A.xyz += [i,j,k]
                    A.xyz /= [l, m, n]
                    A.name = "%s_%i" % (atom.name, dup)
                    Snew.append(A)
                    dup += 1
                # for k
            # for j
        # for i
    # for atom

    return Snew
################################################################################
    
def _testNi():
    a = 3.53
    L = Lattice(a=a,b=a,c=a,alpha=90,beta=90,gamma=90,spacegroup="Fm-3m")
    S = Structure(lattice = L)
    S.title = "Ni"
    A0 = Atom(xyz = numpy.array([0,0,0]), element = "Ni", name = "Ni")
    S.append(A0)
    print S

    # Expand space group
    clist = ['']*10
    clist[0:3] = ["@1", "@2", "@3"]
    Sexp, clistsexp = expandStructureConstraints(S, [clist])
    print "\n****Expanded Cell****"
    print Sexp
    for l in clistsexp: print l

    # Make supercell
    Ssup, clistssup = makeSupercellConstraints(Sexp, clistsexp, 2, 1, 1)
    print "\n****Supercell****"
    print Ssup
    for l in clistssup: print l
    return

def _testGaAs():
    a = 5.6537
    L = Lattice(a=a,b=a,c=a,alpha=90,beta=90,gamma=90,spacegroup="F-43m")
    S = Structure(lattice = L)
    S.title = "GaAs"

    # Expand space group
    clists = []
    A0 = Atom(xyz = numpy.array([0,0,0]), element = "Ga", name = "Ga")
    clist = ["@1", "@2", "@3"]
    clist.extend(['']*7)
    clists.append(clist)
    A1 = Atom(xyz = numpy.array([0.5,0,0]), element = "As", name = "As")
    clist = ["@4", "@5", "@6"]
    clist.extend(['']*7)
    clists.append(clist)
    S.append(A0)
    S.append(A1)
    print S
    Sexp, clistsexp = expandStructureConstraints(S, clists)
    #Sexp.write("temp.stru", "pdffit")
    print "\n****Expanded Cell****"
    print Sexp
    for l in clistsexp: print l

    # Make supercell
    Ssup, clistssup = makeSupercellConstraints(Sexp,clistsexp, 2, 1, 1)
    print "\n****Supercell****"
    print Ssup
    for l in clistssup: print l
    #Ssup.write("temp.stru", "pdffit")
    return

if __name__ == '__main__':
    _testNi()
    _testGaAs()

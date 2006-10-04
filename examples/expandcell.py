#!/usr/bin/env python
from Structure import Structure, Atom, Lattice
from numpy import array

def printSymmetryOperations(S):
    spacegroup = S.lattice.spacegroup
    symops = []
    for symop in\
        spacegroup.symop_list:
            if str(symop) not in symops:
                symops.append(str(symop))
    for symop in symops:
        print symop
    return

def printStructure(S):
    print S.title
    print "space group = %s" % S.lattice.spacegroup.short_name
    print len(S), "atoms"
    for atom in S:
        x, y, z = atom.xyz
        print "%s: %-2.4f %-2.4f %-2.4f" % (atom.name, x, y, z)
    return

def testNi():
    L = Lattice(spacegroup="Fm-3m")
    S = Structure(lattice = L)
    S.title = "Ni"
    A0 = Atom(xyz = array([0,0,0]), element = "Ni", name = "Ni")
    S.append(A0)
    printStructure(S)
    Sexp = S.expandStructure()
    printStructure(Sexp)
    return

def testGaAs():
    L = Lattice(spacegroup="F-43m")
    S = Structure(lattice = L)
    S.title = "GaAs"
    A0 = Atom(xyz = array([0,0,0]), element = "Ga", name = "Ga")
    A1 = Atom(xyz = array([0.25,0.25,0.25]), element = "As", name = "As")
    S.append(A0)
    S.append(A1)
    printStructure(S)
    Sexp = S.expandStructure()
    printStructure(Sexp)
    return (S, Sexp)

if __name__ == "__main__":

    S, Sexp = testGaAs()
    testNi()

    #printSymmetryOperations(S)

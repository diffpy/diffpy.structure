from diffpy.Structure import Structure, StructureFormatError
from diffpy.Structure import Lattice
from diffpy.Structure import Atom

stru = Structure( [ Atom('C', [0,0,0]), Atom('C', [2,2,2]) ],
                lattice=Lattice(3, 3, 3, 90, 90, 90) )
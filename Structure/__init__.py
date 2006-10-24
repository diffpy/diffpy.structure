"""classes related to structure of materials
Classes:
    Atom
    Lattice
    Structure
    PDFFitStructure
Exceptions:
    InvalidStructureFormat
    InvalidLattice
"""

__id__ = "$Id$"

##############################################################################
# interface definitions
##############################################################################

from Structure import Structure, InvalidStructureFormat
from Lattice import Lattice, InvalidLattice
from Atom import Atom
from PDFFitStructure import PDFFitStructure

##############################################################################
# version info
##############################################################################

def _get_module_ids():
    """return list of __id__ values from all modules in this package"""
    import Atom
    import Lattice
    import Structure
    import PDFFitStructure
    import PeriodicTable
    import Parsers
    ids = [
        __id__,
        Atom.__id__,
        Lattice.__id__,
        Structure.__id__,
        PDFFitStructure.__id__,
        PeriodicTable.__id__,
        Parsers._get_package_id()
    ]
    return ids

def _get_package_id():
    """build cumulative ID of this package from module IDs"""
    name = __name__
    id_words = [i.split() for i in _get_module_ids()]
    # are all ids expanded?
    if [ True for idw in id_words if len(idw) < 5 ] != [] :
        package_id = __id__
    # do we have CVS-style ids?
    elif len([ 1 for idw in id_words if "." in idw[2] ]) == len(id_words):
        id_words.sort( lambda x,y : cmp(x[3:5], y[3:5]) )
        date_time_auth = id_words[-1][3:]
        major = id_words[-1][2].split('.')[0]
        minor = sum([int(idw[2].split('.')[1]) for idw in id_words])
        version = "%i.%i" % (major, minor)
        package_id = " ".join(['$Id:', name, version] + date_time_auth)
    # otherwise assume subversion style ids
    else:
        id_words.sort( lambda x,y : cmp(x[2:], y[2:]) )
        last = id_words[-1]
        package_id = " ".join(last[:1]+[name]+last[2:])
    return package_id

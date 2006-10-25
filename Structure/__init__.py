"""classes related to structure of materials
Classes:
    Atom
    Lattice
    Structure
    PDFFitStructure
Exceptions:
    InvalidStructureFormat
    InvalidLattice
    SymmetryError
"""

__id__ = "$Id$"

##############################################################################
# interface definitions
##############################################################################

from Structure import Structure
from Lattice import Lattice
from Atom import Atom
from PDFFitStructure import PDFFitStructure
from exceptions import InvalidStructureFormat, InvalidLattice, SymmetryError

from PeriodicTable import __id__ as aaaa

##############################################################################
# version info
##############################################################################

def _get_module_ids():
    """return list of __id__ values from all modules in this package"""
    import os.path
    import glob
    myFile = os.path.abspath(__file__)
    myDir = os.path.dirname(myFile)
    sourcefiles = glob.glob(os.path.join(myDir, '*.py'))
    modnames = [os.path.splitext(os.path.basename(f))[0] for f in sourcefiles]
    modnames.remove('__init__')
    ids = [__id__]
    for modname in modnames:
        try:
            exec "from %s import __id__ as i" % modname
            ids.append(i)
        except ImportError:
            pass
    return ids

def _get_package_id():
    """build cumulative ID of this package from module IDs"""
    name = __name__
    id_words = [i.split() for i in _get_module_ids()]
    # keep only expanded ids:
    id_words = [idw for idw in id_words if len(idw) >= 5]
    if not id_words:
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

# End of file

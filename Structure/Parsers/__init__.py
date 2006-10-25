"""conversion plugin for various structure formats

The recognized structure formats are defined in submodules, which must be
called P_<format>.  Each of these modules defines Parser class that is
inherited from the StructureParser superclass.  The specific module has
to override parseLines() and toLines() methods of StructureParser.
For a normal usage it should be sufficient to use the parse() and tostring()
methods of this module.
"""

__id__ = "$Id$"

from Structure.Structure import Structure
from Structure.exceptions import InvalidStructureFormat

def _findParsers():
    """return dictionary of recognized formats with associated modules"""
    import sys, os.path
    thisModulePath = os.path.dirname(sys.modules[__name__].__file__)
    from glob import glob
    parserFiles = glob(os.path.join(thisModulePath, 'P_*.py'))
    modules = {}
    for pf in parserFiles:
        pmod = os.path.splitext(os.path.basename(pf))[0]
        pfmt = pmod[2:]
        modules[pfmt] = pmod
    return modules

modules = _findParsers()
formats = modules.keys()
formats.sort()

def parse(s, format):
    """create Structure instance from a string in given format"""
    try:
        pmod = modules[format]
    except KeyError:
        raise InvalidStructureFormat, "no parser for '%s' format" % format
    exec('from ' + pmod + ' import Parser')
    p = Parser()
    lines = s.split('\n')
    if len(lines) and lines[-1] == '':
        del lines[-1]
    stru = p.parseLines(lines)
    return stru

def tostring(stru, format):
    """convert Structure stru to string using the specified format"""
    if not isinstance(stru, Structure):
        raise RuntimeError, "expected instance of Structure"
    try:
        pmod = modules[format]
    except KeyError:
        raise InvalidStructureFormat, "no parser for '%s' format" % format
    exec('from ' + pmod + ' import Parser')
    p = Parser()
    lines = p.toLines(stru)
    s = "\n".join(lines) + "\n"
    return s

def inputFormats():
    """return list of implemented input structure formats"""
    from Structure.Structure import InvalidStructureFormat
    input_formats = []
    for format, module in modules.items():
        try:
            parse("", format)
            input_formats.append(format)
        except InvalidStructureFormat:
            input_formats.append(format)
        except NotImplementedError:
            pass
    input_formats.sort()
    return input_formats

def outputFormats():
    """return list of implemented output structure formats"""
    from Structure.Structure import InvalidStructureFormat
    output_formats = []
    stru = Structure()
    for format, module in modules.items():
        try:
            tostring(stru, format)
            output_formats.append(format)
        except InvalidStructureFormat:
            output_formats.append(format)
        except NotImplementedError:
            pass
    output_formats.sort()
    return output_formats

##############################################################################
# version info
##############################################################################

def _get_module_ids():
    """return list of __id__ values from all modules in this package"""
    import StructureParser
    ids = [ __id__,
            StructureParser.__id__ ]
    for mod in modules.values():
        exec("import %s" % mod)
        ids.append( eval(mod+".__id__") )
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

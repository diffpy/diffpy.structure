########################################################################
#
# Structure         by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2006 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
########################################################################

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

# End of file

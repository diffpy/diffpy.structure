#!/usr/bin/env python

# This script is usually called from top level diffpy installation.

"""Structure - objects for storage and manipulation of crystal structures.

Packages:   diffpy.Structure
Scripts:    transtry, anyeye
"""

# version
__id__ = "$Id$"

from distutils.core import setup
import sys
import os.path

thisfile = os.path.abspath(locals().get('__file__', 'setup.py'))
setup_dir = os.path.dirname(thisfile)

sys.path.insert(0, setup_dir)
from diffpy.Structure import __version__
sys.path.pop(0)
package_version = __version__

# define distribution
setup_args = {
    "name" : "Structure",
    "description" : "Crystal structure container.",
    "version" : package_version,
    "packages" : [ "diffpy.Structure", "diffpy.Structure.Parsers" ],
    "package_dir" : {
        "diffpy.Structure" : os.path.join(setup_dir, "Structure")
        },
    "scripts" : [
        os.path.join(setup_dir, "applications/anyeye"),
        os.path.join(setup_dir, "applications/transtru")
        ],
}

diffpy__init__code = """
import sys
import os.path
mydir = os.path.dirname(__file__)
if mydir not in sys.path:   sys.path.insert(0, mydir)
""".lstrip()

def check_diffpy__init__(distribution):
    """check if diffpy/__init__.py exists and create one if not
    """
    from distutils import log
    if distribution.dry_run:    return
    if 'install_lib' not in distribution.command_obj:   return
    lib_install_dir = distribution.get_command_obj('install_lib').install_dir
    initfile = os.path.join(lib_install_dir, 'diffpy', '__init__.py')
    if os.path.isfile(initfile):    return
    # we need to create and compile the file
    log.info("creating " + initfile)
    out = open(initfile, 'w')
    out.write(diffpy__init__code)
    out.close()
    import compiler
    log.info("byte-compiling %s to %s" % \
            (initfile, os.path.basename(initfile)) )
    compiler.compileFile(initfile)
    return

if __name__ == "__main__":
    distribution = setup(**setup_args)
    check_diffpy__init__(distribution)

# End of file 

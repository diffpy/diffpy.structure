#!/usr/bin/env python

# version
__id__ = "$Id$"

from distutils.core import setup
import sys
import os.path

thisfile = os.path.abspath(locals().get('__file__', 'setup.py'))
setup_dir = os.path.dirname(thisfile)

sys.path.insert(0, setup_dir)
from Structure import __version__
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
    """check if diffpy has __init__.py and create one if not
    """
    from distutils import log
    install_lib = None
    if 'install' in distribution.commands:
        opts = distribution.get_option_dict('install')
        install_lib = opts.get('install_lib', 2*[None])[1]
    if 'install_lib' in distribution.commands and not install_lib:
        opts = distribution.get_option_dict('install_lib')
        install_lib = opts.get('install_dir', 2*[None])[1]
    if not install_lib:             return
    initfile = os.path.join(install_lib, 'diffpy', '__init__.py')
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

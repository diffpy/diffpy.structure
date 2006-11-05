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

if __name__ == "__main__":
    distribution = setup(**setup_args)

# End of file 

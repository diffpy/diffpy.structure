# This module is imported from top level diffpy setup.py.
# It has to define the following variables:
#     name, description, diffpy_deps, other_deps, setup_args
# Optional variables:
#     makefiles -- a list of Makefiles to be build before installation

"""Structure - objects for storage and manipulation of crystal structures.

Packages:   diffpy.Structure
Scripts:    transtry, anyeye
"""

# version
__id__ = "$Id$"

import os.path

thisfile = os.path.abspath(locals().get('__file__', 'setup_args.py'))
thisdir = os.path.dirname(thisfile)

def prependThisDir(files):
    return [os.path.join(thisdir, f) for f in files]

# name of this subpackage
name = "diffpy.Structure"
description =  "Crystal structure container.",

# dependencies from diffpy
diffpy_deps = []

# third-party dependencies
other_deps = [
    "numpy",
    "CifFile",
    ]

# define distribution arguments for this subpackage
setup_args = {
    "name" : name,
    "description" : description,
    "packages" : [
        "diffpy.Structure",
        "diffpy.Structure.Parsers",
        "diffpy.Structure.expansion",
        ],
    "package_dir" : {
        "diffpy.Structure" : os.path.join(thisdir, "Structure"),
        },
#   "scripts" : prependThisDir([
#       "applications/anyeye",
#       "applications/transtru",
#       ]),
}

# End of file 

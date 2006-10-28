#!/usr/bin/env python

from distutils.core import setup
import sys

import Structure.version
package_version = Structure.version

# define distribution
dist = setup(
    name = "Structure",
    description = "Crystal structure container.",
    version = package_version,
    packages = [ "diffpy.Structure", "diffpy.Structure.Parsers" ],
    package_dir = { "diffpy.Structure" : "Structure" },
    scripts = [ "applications/anyeye", "applications/transtru" ]
)

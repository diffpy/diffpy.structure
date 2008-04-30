#!/usr/bin/env python

# Installation script for diffpy.Structure

"""Structure - objects for storage and manipulation of crystal structures.

Packages:   diffpy.Structure
Scripts:    transtry, anyeye
"""

# version
__id__ = "$Id$"

from setuptools import setup, find_packages
import fix_setuptools_chmod

# define distribution
dist = setup(
        name = "diffpy.Structure",
        namespace_packages = ['diffpy'],
        version = "1.0c1",
        packages = [
            "diffpy",
            "diffpy.Structure",
            "diffpy.Structure.Parsers",
            "diffpy.Structure.expansion",
        ],
        install_requires = ['PyCifRW'],
        dependency_links = [
            "http://diffpy.org/packages/",
        ],

        author = "Simon J.L. Billinge",
        author_email = "sb2896@columbia.edu",
        description = "Empty namespace module for the diffpy library.",
        license = "BSD",
        url = "http://www.diffpy.org/",
        keywords = "diffpy Structure container",
)

# End of file

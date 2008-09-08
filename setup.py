#!/usr/bin/env python

# Installation script for diffpy.Structure

"""Structure - objects for storage and manipulation of crystal structures.

Packages:   diffpy.Structure
Scripts:    transtry, anyeye
"""

from setuptools import setup, find_packages
import fix_setuptools_chmod

# define distribution
dist = setup(
        name = "diffpy.Structure",
        version = "1.0c1",
        namespace_packages = ['diffpy'],
        packages = find_packages(),
        install_requires = ['PyCifRW'],
        dependency_links = [
            "http://www.diffpy.org/packages/",
        ],

        author = "Simon J.L. Billinge",
        author_email = "sb2896@columbia.edu",
        description = "Crystal structure container " + \
                      "and parsers for structure formats.",
        license = "BSD",
        url = "http://www.diffpy.org/",
        keywords = "diffpy Structure container",
)

# End of file

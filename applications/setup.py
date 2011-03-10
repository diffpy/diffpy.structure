#!/usr/bin/env python

# This script installs entry points to diffpy.Structure.applications
# It does not install any Python packages.

"""Scripts:

anyeye   -- wrapper for structure viewer that supports extra file formats.
transtru -- translator of structure formats supported by diffpy.Structure
"""

from setuptools import setup

# define distribution
setup(
        name = 'diffpy.Structure.applications',
        version = '1.1',
        entry_points = {
            'console_scripts': [
                'anyeye=diffpy.Structure.applications.anyeye:main',
                'transtru=diffpy.Structure.applications.transtru:main',
            ],
        },
        install_requires = [
            'diffpy.Structure>=1.1',
        ],
        dependency_links = [
            # FIXME: remove dev.danse.us for final release
            'http://dev.danse.us/packages/',
            'http://www.diffpy.org/packages/',
        ],

        author = 'Simon J.L. Billinge',
        author_email = 'sb2896@columbia.edu',
        maintainer = 'Pavol Juhas',
        maintainer_email = 'pj2192@columbia.edu',
        url = 'http://www.diffpy.org/',
        download_url = 'http://www.diffpy.org/packages/',
        description =
            "Structure viewer wrapper that supports more file formats.",
        license = 'BSD',
        keywords = 'structure viewer format support',
        classifiers = [
            # List of possible values at
            # http://pypi.python.org/pypi?:action=list_classifiers
            'Development Status :: 5 - Production/Stable',
            'Environment :: MacOS X',
            'Environment :: Win32 (MS Windows)',
            'Intended Audience :: Science/Research',
            'Operating System :: MacOS',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: POSIX',
            'Programming Language :: Python :: 2.5',
            'Topic :: Scientific/Engineering :: Physics',
        ],
)

# End of file

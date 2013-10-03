#!/usr/bin/env python

# Installation script for diffpy.Structure

"""Structure - objects for storage and manipulation of crystal structures.

Packages:   diffpy.Structure
"""

import os
from setuptools import setup, find_packages

def gitversion():
    from subprocess import Popen, PIPE
    proc = Popen(['git', 'describe'], stdout=PIPE)
    desc = proc.stdout.read().strip()
    proc = Popen(['git', 'log', '-1', '--format=%ai'], stdout=PIPE)
    isodate = proc.stdout.read()
    date = isodate.split()[0].replace('-', '')
    rv = desc + '-' + date
    return rv


def getsetupcfg():
    cfgfile = 'setup.cfg'
    from ConfigParser import SafeConfigParser
    cp = SafeConfigParser()
    cp.read(cfgfile)
    if not os.path.isdir('.git'):  return cp
    d = cp.defaults()
    vcfg = d.get('version', '')
    vgit = gitversion()
    if vgit != vcfg:
        cp.set('DEFAULT', 'version', vgit)
        cp.write(open(cfgfile, 'w'))
    return cp

cp = getsetupcfg()

# define distribution
setup(
        name = "diffpy.Structure",
        version = cp.get('DEFAULT', 'version'),
        namespace_packages = ['diffpy'],
        packages = find_packages(),
        test_suite = 'diffpy.Structure.tests',
        include_package_data = True,
        install_requires = [
            'PyCifRW',
        ],
        dependency_links = [
            'http://www.diffpy.org/packages/',
        ],

        author = 'Simon J.L. Billinge',
        author_email = 'sb2896@columbia.edu',
        maintainer = 'Pavol Juhas',
        maintainer_email = 'pj2192@columbia.edu',
        url = 'http://www.diffpy.org/',
        download_url = 'http://www.diffpy.org/packages/',
        description = "Crystal structure container " + \
                      "and parsers for structure formats.",
        license = 'BSD',
        keywords = "diffpy Structure container",
        classifiers = [
            # List of possible values at
            # http://pypi.python.org/pypi?:action=list_classifiers
            'Development Status :: 5 - Production/Stable',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'Operating System :: MacOS',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: POSIX',
            'Programming Language :: Python :: 2.5',
            'Topic :: Scientific/Engineering :: Physics',
        ],
)

# End of file

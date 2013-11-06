#!/usr/bin/env python

# Installation script for diffpy.Structure

"""Structure - objects for storage and manipulation of crystal structures.

Packages:   diffpy.Structure
"""

import os
from setuptools import setup, find_packages

# versioncfgfile holds version data for git commit hash and date.
# It must reside in the same directory as version.py.
versioncfgfile = 'diffpy/Structure/version.cfg'

def gitinfo():
    from subprocess import Popen, PIPE
    proc = Popen(['git', 'describe'], stdout=PIPE)
    desc = proc.stdout.read()
    proc = Popen(['git', 'log', '-1', '--format=%H %ai'], stdout=PIPE)
    glog = proc.stdout.read()
    rv = {}
    rv['version'] = '-'.join(desc.strip().split('-')[:2])
    rv['commit'], rv['date'] = glog.strip().split(None, 1)
    return rv


def getversioncfg():
    import os
    from ConfigParser import SafeConfigParser
    cp = SafeConfigParser()
    cp.read(versioncfgfile)
    if not os.path.isdir('.git'):  return cp
    d = cp.defaults()
    g = gitinfo()
    if g['commit'] != d.get('commit'):
        cp.set('DEFAULT', 'version', g['version'])
        cp.set('DEFAULT', 'commit', g['commit'])
        cp.set('DEFAULT', 'date', g['date'])
        cp.write(open(versioncfgfile, 'w'))
    return cp

cp = getversioncfg()

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

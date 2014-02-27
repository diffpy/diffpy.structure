#!/usr/bin/env python

# Installation script for diffpy.Structure

"""Structure - objects for storage and manipulation of crystal structures.

Packages:   diffpy.Structure
"""

import os
from setuptools import setup, find_packages

# versioncfgfile holds version data for git commit hash and date.
# It must reside in the same directory as version.py.
MYDIR = os.path.dirname(os.path.abspath(__file__))
versioncfgfile = os.path.join(MYDIR, 'diffpy/Structure/version.cfg')

def gitinfo():
    from subprocess import Popen, PIPE
    kw = dict(stdout=PIPE, cwd=MYDIR)
    proc = Popen(['git', 'describe', '--match=v[[:digit:]]*'], **kw)
    desc = proc.stdout.read()
    proc = Popen(['git', 'log', '-1', '--format=%H %at %ai'], **kw)
    glog = proc.stdout.read()
    rv = {}
    rv['version'] = '-'.join(desc.strip().split('-')[:2]).lstrip('v')
    rv['commit'], rv['timestamp'], rv['date'] = glog.strip().split(None, 2)
    return rv


def getversioncfg():
    from ConfigParser import SafeConfigParser
    cp = SafeConfigParser()
    cp.read(versioncfgfile)
    gitdir = os.path.join(MYDIR, '.git')
    if not os.path.isdir(gitdir):  return cp
    d = cp.defaults()
    g = gitinfo()
    if g['version'] != d.get('version') or g['commit'] != d.get('commit'):
        cp.set('DEFAULT', 'version', g['version'])
        cp.set('DEFAULT', 'commit', g['commit'])
        cp.set('DEFAULT', 'date', g['date'])
        cp.set('DEFAULT', 'timestamp', g['timestamp'])
        cp.write(open(versioncfgfile, 'w'))
    return cp

versiondata = getversioncfg()

# define distribution
setup_args = dict(
        name = "diffpy.Structure",
        version = versiondata.get('DEFAULT', 'version'),
        namespace_packages = ['diffpy'],
        packages = find_packages(),
        test_suite = 'diffpy.Structure.tests',
        include_package_data = True,
        zip_safe = False,
        install_requires = [
            'PyCifRW',
        ],

        author = 'Simon J.L. Billinge',
        author_email = 'sb2896@columbia.edu',
        maintainer = 'Pavol Juhas',
        maintainer_email = 'pavol.juhas@gmail.com',
        url = 'https://github.com/diffpy/diffpy.Structure',
        #download_url = 'http://www.diffpy.org/packages/',
        description = "Crystal structure container " + \
                      "and parsers for structure formats.",
        license = 'BSD, see LICENSE.txt',
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
            'Programming Language :: Python :: 2.6',
            'Programming Language :: Python :: 2.7'
            'Topic :: Scientific/Engineering :: Chemistry',
            'Topic :: Scientific/Engineering :: Physics',
        ],
)

if __name__ == '__main__':
    setup(**setup_args)

# End of file

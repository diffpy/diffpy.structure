#!/usr/bin/env python

"""Objects for storage and manipulation of crystal structure data.

Packages:   diffpy.Structure
"""

import os
import re
import sys
from setuptools import setup, find_packages

# Use this version when git data are not available, like in git zip archive.
# Update when tagging a new release.
FALLBACK_VERSION = '1.3.5.post0'

# versioncfgfile holds version data for git commit hash and date.
# It must reside in the same directory as version.py.
MYDIR = os.path.dirname(os.path.abspath(__file__))
versioncfgfile = os.path.join(MYDIR, 'src/diffpy/Structure/version.cfg')
gitarchivecfgfile = versioncfgfile.replace('version.cfg', 'gitarchive.cfg')


def gitinfo():
    from subprocess import Popen, PIPE
    kw = dict(stdout=PIPE, cwd=MYDIR, universal_newlines=True)
    proc = Popen(['git', 'describe', '--match=v[[:digit:]]*'], **kw)
    desc = proc.stdout.read()
    proc = Popen(['git', 'log', '-1', '--format=%H %at %ai'], **kw)
    glog = proc.stdout.read()
    rv = {}
    rv['version'] = '.post'.join(desc.strip().split('-')[:2]).lstrip('v')
    rv['commit'], rv['timestamp'], rv['date'] = glog.strip().split(None, 2)
    return rv


def getversioncfg():
    if sys.version_info[0] >= 3:
        from configparser import RawConfigParser
    else:
        from ConfigParser import RawConfigParser
    vd0 = dict(version=FALLBACK_VERSION, commit='', date='', timestamp=0)
    # first fetch data from gitarchivecfgfile, ignore if it is unexpanded
    g = vd0.copy()
    cp0 = RawConfigParser(vd0)
    cp0.read(gitarchivecfgfile)
    if '$Format:' not in cp0.get('DEFAULT', 'commit'):
        g = cp0.defaults()
        mx = re.search(r'\btag: v(\d[^,]*)', g.pop('refnames'))
        if mx:
            g['version'] = mx.group(1)
    # then try to obtain version data from git.
    gitdir = os.path.join(MYDIR, '.git')
    if os.path.exists(gitdir) or 'GIT_DIR' in os.environ:
        try:
            g = gitinfo()
        except OSError:
            pass
    # finally, check and update the active version file
    cp = RawConfigParser()
    cp.read(versioncfgfile)
    d = cp.defaults()
    rewrite = not d or (g['commit'] and (
        g['version'] != d.get('version') or g['commit'] != d.get('commit')))
    if rewrite:
        cp.set('DEFAULT', 'version', g['version'])
        cp.set('DEFAULT', 'commit', g['commit'])
        cp.set('DEFAULT', 'date', g['date'])
        cp.set('DEFAULT', 'timestamp', g['timestamp'])
        with open(versioncfgfile, 'w') as fp:
            cp.write(fp)
    return cp

versiondata = getversioncfg()

# define distribution
setup_args = dict(
    name = "diffpy.Structure",
    version = versiondata.get('DEFAULT', 'version'),
    namespace_packages = ['diffpy'],
    packages = find_packages('src'),
    package_dir = {'' : 'src'},
    test_suite = 'diffpy.Structure.tests',
    include_package_data = True,
    zip_safe = False,
    install_requires = [
        'pycifrw',
    ],

    author = 'Simon J.L. Billinge group',
    author_email = 'sb2896@columbia.edu',
    maintainer = 'Pavol Juhas',
    maintainer_email = 'pavol.juhas@gmail.com',
    url = 'https://github.com/diffpy/diffpy.structure',
    description = "Crystal structure container "
                  "and parsers for structure formats.",
    license = 'BSD-style license',
    keywords = "crystal Structure data storage CIF PDB",
    classifiers = [
        # List of possible values at
        # http://pypi.python.org/pypi?:action=list_classifiers
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
    ],
)

if __name__ == '__main__':
    setup(**setup_args)

# End of file

#!/usr/bin/env python

from distutils.core import setup
from distutils import log
import sys, os, os.path

if "-h" in sys.argv or "--help" in sys.argv:
    print "Special option:"
    print "  --slapper           set options reasonable for slapper"
    print

slapper_root = "/u24/local"
slapper_lib_dir = os.path.join(slapper_root,
        "bgPython/python%i.%i-packages" % sys.version_info[:2] )
slapper_scripts_dir = os.path.join(slapper_root, "bgPython/applications")
slapper_bin_dir = os.path.join(slapper_root, "bin")
slapper_hashbang = "#!/usr/bin/python"

# check if we are installing on slapper
if "--slapper" in sys.argv:
    slapper_installation = True
    sys.argv.remove("--slapper")
    options = { "install" : {
                   "install_lib"     : slapper_lib_dir,
                   "install_scripts" : slapper_scripts_dir }
              }
else:
    slapper_installation = False
    options = {"install" : {}}

# define distribution
dist = setup(
    name = "Structure",
    description = "Objects for storage and manipulations of atom structures",
    packages = [ "Structure", "Structure.Parsers" ],
#   package_dir = { "Structure" : "", "Structure.Parsers" : "Parsers" },
    options = options,
    scripts = [ "applications/anyeye", "applications/transtru" ]
)

# make symbolic links and adjust hashbang if we are installing on slapper
if slapper_installation \
and dist.commands in [ ['install'], ['install_scripts'] ] \
and dist.scripts is not None:
    for f in dist.scripts:
        src = os.path.join('../bgPython/', f)
        dest = os.path.join(slapper_bin_dir, os.path.basename(f))
        if not os.path.isfile(dest):
            log.info( "making symlink %s --> %s" % (dest, src) )
            os.symlink(src, dest)
    for f in dist.scripts:
        f_full = os.path.join(slapper_scripts_dir, os.path.basename(f))
        f_lines = open(f_full).readlines()
        if f_lines[0].strip() != slapper_hashbang:
            f_lines[0] = slapper_hashbang + "\n"
            open(f_full,'w').writelines(f_lines)
            log.info( "adjusted hashbang in " + f_full )

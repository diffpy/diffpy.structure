#!/usr/bin/env python
##############################################################################
#
# diffpy.structure  by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2006 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################

"""Translate structure file to different format.

Usage: ``transtru INFMT..OUTFMT strufile``

Translates structure file strufile from `INFMT` to `OUTFMT` format and prints it
to the screen. Use "-" as `strufile` to read from standard input. To save the
translated file, use

    ``transtru INFMT..OUTFMT strufile > strufile.out``

Supported input and output structure formats are
    * `INFMT`: ``inputFormats``
    * `OUTFMT`: ``outputFormats``

Options:
  -h, --help
    Display this message.

  -V, --version
    Show script version.
"""

from __future__ import print_function

import sys

from diffpy.structure import Structure
from diffpy.structure.structureerrors import StructureFormatError


def usage(style=None):
    """Show usage info, for ``style=="brief"`` show only first 2 lines."""
    import os.path

    myname = os.path.basename(sys.argv[0])
    msg = __doc__.replace("transtru", myname)
    if style == "brief":
        msg = msg.split("\n")[1] + "\n" + "Try `%s --help' for more information." % myname
    else:
        from diffpy.structure.parsers import inputFormats, outputFormats

        msg = msg.replace("inputFormats", " ".join(inputFormats()))
        msg = msg.replace("outputFormats", " ".join(outputFormats()))
    print(msg)
    return


def version():
    from diffpy.structure import __version__

    print("diffpy.structure", __version__)
    return


def main():
    import getopt

    # default parameters
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hV", ["help", "version"])
    except getopt.GetoptError as errmsg:
        print(errmsg, file=sys.stderr)
        sys.exit(2)
    # process options
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-V", "--version"):
            version()
            sys.exit()
    if len(args) < 1:
        usage("brief")
        sys.exit()
    # process arguments
    from diffpy.structure.parsers import inputFormats, outputFormats

    try:
        infmt, outfmt = args[0].split("..", 1)
        if infmt not in inputFormats():
            print("'%s' is not valid input format" % infmt, file=sys.stderr)
            sys.exit(2)
        if outfmt not in outputFormats():
            print("'%s' is not valid output format" % outfmt, file=sys.stderr)
            sys.exit(2)
    except ValueError:
        print("invalid format specification '%s' does not contain .." % args[0], file=sys.stderr)
        sys.exit(2)
    # ready to do some real work
    try:
        strufile = args[1]
        stru = Structure()
        if args[1] == "-":
            stru.readStr(sys.stdin.read(), infmt)
        else:
            stru.read(strufile, infmt)
        sys.stdout.write(stru.writeStr(outfmt))
    except IndexError:
        print("strufile not specified", file=sys.stderr)
        sys.exit(2)
    except IOError as e:
        print("%s: %s" % (strufile, e.strerror), file=sys.stderr)
        sys.exit(1)
    except StructureFormatError as e:
        print("%s: %s" % (strufile, e), file=sys.stderr)
        sys.exit(1)
    return


if __name__ == "__main__":
    main()

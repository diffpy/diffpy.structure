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

"""
Anyeye view structure file in atomeye.

Usage: ``anyeye [options] strufile``

Anyeye understands more `Structure` formats than atomeye. It converts `strufile`
to a temporary XCFG file which is opened in atomeye. See supported file formats:
``inputFormats``

Options:
  -f, --formula
      Override chemical formula in `strufile`. The formula defines
      elements in the same order as in `strufile`, e.g., ``Na4Cl4``.

  -w, --watch
      Watch input file for changes.

  --viewer=VIEWER
      The structure viewer program, by default "atomeye".
      The program will be executed as "VIEWER structurefile".

  --formats=FORMATS
      Comma-separated list of file formats that are understood
      by the VIEWER, by default ``"xcfg,pdb"``. Files of other
      formats will be converted to the first listed format.

  -h, --help
      Display this message and exit.

  -V, --version
      Show script version and exit.
"""

from __future__ import print_function

import os
import re
import signal
import sys

from diffpy.structure.structureerrors import StructureFormatError

# parameter dictionary
pd = {
    "formula": None,
    "watch": False,
    "viewer": "atomeye",
    "formats": ["xcfg", "pdb"],
}


def usage(style=None):
    """Show usage info, for ``style=="brief"`` show only first 2 lines."""
    import os.path

    myname = os.path.basename(sys.argv[0])
    msg = __doc__.replace("anyeye", myname)
    if style == "brief":
        msg = msg.split("\n")[1] + "\n" + "Try `%s --help' for more information." % myname
    else:
        from diffpy.structure.parsers import inputFormats

        fmts = [f for f in inputFormats() if f != "auto"]
        msg = msg.replace("inputFormats", " ".join(fmts))
    print(msg)
    return


def version():
    from diffpy.structure import __version__

    print("anyeye", __version__)
    return


def loadStructureFile(filename, format="auto"):
    """Load structure from specified file.

    Parameters
    ----------
    filename : str
        Path to the structure file.
    format : str, Optional
        File format, by default "auto".

    Returns
    -------
    tuple
        A tuple of (Structure, fileformat).
    """
    from diffpy.structure import Structure

    stru = Structure()
    p = stru.read(filename, format)
    fileformat = p.format
    return (stru, fileformat)


def convertStructureFile(pd):
    # make temporary directory on the first pass
    if "tmpdir" not in pd:
        from tempfile import mkdtemp

        pd["tmpdir"] = mkdtemp()
    strufile = pd["strufile"]
    tmpfile = os.path.join(pd["tmpdir"], os.path.basename(strufile))
    pd["tmpfile"] = tmpfile
    # speed up file processing in the watch mode
    fmt = pd.get("format", "auto")
    stru = None
    if fmt == "auto":
        stru, fmt = loadStructureFile(strufile)
        pd["fmt"] = fmt
    # if fmt is recognized by the viewer, use as is
    if fmt in pd["formats"] and pd["formula"] is None:
        import shutil

        shutil.copyfile(strufile, tmpfile + ".tmp")
        os.rename(tmpfile + ".tmp", tmpfile)
        return
    # otherwise convert to the first recognized viewer format
    if stru is None:
        stru = loadStructureFile(strufile, fmt)[0]
    if pd["formula"]:
        formula = pd["formula"]
        if len(formula) != len(stru):
            emsg = "Formula has %i atoms while structure %i" % (len(formula), len(stru))
            raise RuntimeError(emsg)
        for a, el in zip(stru, formula):
            a.element = el
    elif format == "rawxyz":
        for a in stru:
            if a.element == "":
                a.element = "C"
    stru.write(tmpfile + ".tmp", pd["formats"][0])
    os.rename(tmpfile + ".tmp", tmpfile)
    return


def watchStructureFile(pd):
    from time import sleep

    strufile = pd["strufile"]
    tmpfile = pd["tmpfile"]
    while pd["watch"]:
        if os.path.getmtime(tmpfile) < os.path.getmtime(strufile):
            convertStructureFile(pd)
        sleep(1)
    return


def cleanUp(pd):
    if "tmpfile" in pd:
        os.remove(pd["tmpfile"])
        del pd["tmpfile"]
    if "tmpdir" in pd:
        os.rmdir(pd["tmpdir"])
        del pd["tmpdir"]
    return


def parseFormula(formula):
    """Parse chemical formula and return a list of elements"""
    # remove all blanks
    formula = re.sub(r"\s", "", formula)
    if not re.match("^[A-Z]", formula):
        raise RuntimeError("InvalidFormula '%s'" % formula)
    elcnt = re.split("([A-Z][a-z]?)", formula)[1:]
    ellst = []
    try:
        for i in range(0, len(elcnt), 2):
            el = elcnt[i]
            cnt = elcnt[i + 1]
            cnt = (cnt == "") and 1 or int(cnt)
            ellst.extend(cnt * [el])
    except ValueError:
        emsg = "Invalid formula, %r is not valid count" % elcnt[i + 1]
        raise RuntimeError(emsg)
    return ellst


def die(exit_status=0, pd={}):
    cleanUp(pd)
    sys.exit(exit_status)


def signalHandler(signum, stackframe):
    # revert to default handler
    signal.signal(signum, signal.SIG_DFL)
    if signum == signal.SIGCHLD:
        pid, exit_status = os.wait()
        exit_status = (exit_status >> 8) + (exit_status & 0x00FF)
        die(exit_status, pd)
    else:
        die(1, pd)
    return


def main():
    import getopt

    # default parameters
    pd["watch"] = False
    try:
        opts, args = getopt.getopt(
            sys.argv[1:], "f:whV", ["formula=", "watch", "viewer=", "formats=", "help", "version"]
        )
    except getopt.GetoptError as errmsg:
        print(errmsg, file=sys.stderr)
        die(2)
    # process options
    for o, a in opts:
        if o in ("-f", "--formula"):
            try:
                pd["formula"] = parseFormula(a)
            except RuntimeError as msg:
                print(msg, file=sys.stderr)
                die(2)
        elif o in ("-w", "--watch"):
            pd["watch"] = True
        elif o == "--viewer":
            pd["viewer"] = a
        elif o == "--formats":
            pd["formats"] = [w.strip() for w in a.split(",")]
        elif o in ("-h", "--help"):
            usage()
            die()
        elif o in ("-V", "--version"):
            version()
            die()
    if len(args) < 1:
        usage("brief")
        die()
    elif len(args) > 1:
        print("too many structure files", file=sys.stderr)
        die(2)
    pd["strufile"] = args[0]
    # trap the following signals
    signal.signal(signal.SIGHUP, signalHandler)
    signal.signal(signal.SIGQUIT, signalHandler)
    signal.signal(signal.SIGSEGV, signalHandler)
    signal.signal(signal.SIGTERM, signalHandler)
    signal.signal(signal.SIGINT, signalHandler)
    env = os.environ.copy()
    if os.path.basename(pd["viewer"]).startswith("atomeye"):
        env["XLIB_SKIP_ARGB_VISUALS"] = "1"
    # try to run the thing:
    try:
        convertStructureFile(pd)
        spawnargs = (pd["viewer"], pd["viewer"], pd["tmpfile"], env)
        # load strufile in atomeye
        if pd["watch"]:
            signal.signal(signal.SIGCLD, signalHandler)
            os.spawnlpe(os.P_NOWAIT, *spawnargs)
            watchStructureFile(pd)
        else:
            status = os.spawnlpe(os.P_WAIT, *spawnargs)
            die(status, pd)
    except IOError as e:
        print("%s: %s" % (args[0], e.strerror), file=sys.stderr)
        die(1, pd)
    except StructureFormatError as e:
        print("%s: %s" % (args[0], e), file=sys.stderr)
        die(1, pd)
    return


if __name__ == "__main__":
    main()

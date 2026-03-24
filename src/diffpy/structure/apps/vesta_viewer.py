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
"""View structure file in VESTA.

Usage: ``vestaview [options] strufile``

Vestaview understands more `Structure` formats than VESTA. It converts
`strufile` to a temporary VESTA or CIF file which is opened in VESTA.
See supported file formats: ``inputFormats``

Options:
  -f, --formula
      Override chemical formula in `strufile`. The formula defines
      elements in the same order as in `strufile`, e.g., ``Na4Cl4``.

  -w, --watch
      Watch input file for changes.

  --viewer=VIEWER
      The structure viewer program, by default "vesta".
      The program will be executed as "VIEWER structurefile".

  --formats=FORMATS
      Comma-separated list of file formats that are understood
      by the VIEWER, by default ``"vesta,cif"``. Files of other
      formats will be converted to the first listed format.

  -h, --help
      Display this message and exit.

  -V, --version
      Show script version and exit.

Notes
-----
VESTA is the actively maintained successor to AtomEye. Unlike AtomEye,
VESTA natively reads CIF, its own ``.vesta`` format, and several other
crystallographic file types, so format conversion is only required for
formats not in that set.

AtomEye XCFG format is no longer a default target format but the XCFG
parser (``P_xcfg``) remains available in ``diffpy.structure.parsers``
for backward compatibility.
"""

from __future__ import print_function

import os
import re
import signal
import sys

from diffpy.structure.structureerrors import StructureFormatError

pd = {
    "formula": None,
    "watch": False,
    "viewer": "vesta",
    "formats": ["vesta", "cif"],
}


def usage(style=None):
    """Show usage info; for ``style=="brief"`` show only first 2
    lines."""
    import os.path

    myname = os.path.basename(sys.argv[0])
    msg = __doc__.replace("vestaview", myname)
    if style == "brief":
        msg = msg.split("\n")[1] + "\n" + "Try `%s --help' for more information." % myname
    else:
        from diffpy.structure.parsers import input_formats

        fmts = [f for f in input_formats() if f != "auto"]
        msg = msg.replace("inputFormats", " ".join(fmts))
    print(msg)
    return


def version():
    from diffpy.structure import __version__

    print("vestaview", __version__)
    return


def load_structure_file(filename, format="auto"):
    """Load structure from specified file.

    Parameters
    ----------
    filename : str
        Path to the structure file.
    format : str, optional
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


def convert_structure_file(pd):
    """Convert `strufile` to a temporary file understood by the viewer.

    On first call a temporary directory is created and stored in *pd*.
    Subsequent calls in watch mode reuse the directory.

    The VESTA viewer natively reads ``.vesta`` and ``.cif`` files, so if
    the source is already in one of the formats listed in ``pd["formats"]``
    and no formula override is requested the file is copied unchanged.
    Otherwise the structure is loaded and re-written in the first format
    listed in ``pd["formats"]``.

    Parameters
    ----------
    pd : dict
        Parameter dictionary containing at minimum ``"strufile"`` and
        ``"formats"`` keys.  Modified in-place to add ``"tmpdir"`` and
        ``"tmpfile"`` on first call.
    """
    # Make temporary directory on the first pass.
    if "tmpdir" not in pd:
        from tempfile import mkdtemp

        pd["tmpdir"] = mkdtemp()
    strufile = pd["strufile"]
    tmpfile = os.path.join(pd["tmpdir"], os.path.basename(strufile))
    pd["tmpfile"] = tmpfile
    # Speed up file processing in the watch mode by caching format.
    fmt = pd.get("format", "auto")
    stru = None
    if fmt == "auto":
        stru, fmt = load_structure_file(strufile)
        pd["fmt"] = fmt
    # If fmt is already recognised by the viewer and no override, copy as-is.
    if fmt in pd["formats"] and pd["formula"] is None:
        import shutil

        shutil.copyfile(strufile, tmpfile + ".tmp")
        os.rename(tmpfile + ".tmp", tmpfile)
        return
    # Otherwise convert to the first viewer-recognised format.
    if stru is None:
        stru = load_structure_file(strufile, fmt)[0]
    if pd["formula"]:
        formula = pd["formula"]
        if len(formula) != len(stru):
            emsg = "Formula has %i atoms while structure %i" % (
                len(formula),
                len(stru),
            )
            raise RuntimeError(emsg)
        for a, el in zip(stru, formula):
            a.element = el
    elif fmt == "rawxyz":
        for a in stru:
            if a.element == "":
                a.element = "C"
    stru.write(tmpfile + ".tmp", pd["formats"][0])
    os.rename(tmpfile + ".tmp", tmpfile)
    return


def watch_structure_file(pd):
    """Watch *strufile* for modifications and reconvert when changed.

    Polls the modification timestamps of ``pd["strufile"]`` and
    ``pd["tmpfile"]`` once per second. When the source is newer the
    file is reconverted via :func:`convert_structure_file`.

    Parameters
    ----------
    pd : dict
        Parameter dictionary as used by :func:`convert_structure_file`.
    """
    from time import sleep

    strufile = pd["strufile"]
    tmpfile = pd["tmpfile"]
    while pd["watch"]:
        if os.path.getmtime(tmpfile) < os.path.getmtime(strufile):
            convert_structure_file(pd)
        sleep(1)
    return


def clean_up(pd):
    """Remove temporary file and directory created by
    :func:`convert_structure_file`.

    Parameters
    ----------
    pd : dict
        Parameter dictionary that may contain ``"tmpfile"`` and
        ``"tmpdir"`` entries to be removed.
    """
    if "tmpfile" in pd:
        os.remove(pd["tmpfile"])
        del pd["tmpfile"]
    if "tmpdir" in pd:
        os.rmdir(pd["tmpdir"])
        del pd["tmpdir"]
    return


def parse_formula(formula):
    """Parse chemical formula and return a list of elements.

    Parameters
    ----------
    formula : str
        Chemical formula string such as ``"Na4Cl4"`` or ``"H2O"``.

    Returns
    -------
    list of str
        Ordered list of element symbols with repetition matching the
        formula, e.g. ``["Na", "Na", "Na", "Na", "Cl", "Cl", "Cl", "Cl"]``.

    Raises
    ------
    RuntimeError
        When *formula* does not start with an uppercase letter or contains
        a non-integer count.
    """
    # Remove all whitespace.
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
    """Clean up temporary files and exit with *exit_status*.

    Parameters
    ----------
    exit_status : int, optional
        Exit code passed to :func:`sys.exit`, by default 0.
    pd : dict, optional
        Parameter dictionary forwarded to :func:`clean_up`.
    """
    clean_up(pd)
    sys.exit(exit_status)


def signal_handler(signum, stackframe):
    """Handle OS signals by reverting to the default handler and
    exiting.

    On ``SIGCHLD`` the child exit status is harvested via
    :func:`os.wait`; on all other signals :func:`die` is called with
    exit status 1.

    Parameters
    ----------
    signum : int
        Signal number.
    stackframe : frame
        Current stack frame (unused).
    """
    # Revert to default handler before acting to avoid re-entrancy.
    signal.signal(signum, signal.SIG_DFL)
    if signum == signal.SIGCHLD:
        pid, exit_status = os.wait()
        exit_status = (exit_status >> 8) + (exit_status & 0x00FF)
        die(exit_status, pd)
    else:
        die(1, pd)
    return


def main():
    """Entry point for the ``vestaview`` command-line tool."""
    import getopt

    # Reset to defaults each invocation.
    pd["watch"] = False
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "f:whV",
            ["formula=", "watch", "viewer=", "formats=", "help", "version"],
        )
    except getopt.GetoptError as errmsg:
        print(errmsg, file=sys.stderr)
        die(2)
    # Process options.
    for o, a in opts:
        if o in ("-f", "--formula"):
            try:
                pd["formula"] = parse_formula(a)
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
    # Trap the following signals.
    signal.signal(signal.SIGHUP, signal_handler)
    signal.signal(signal.SIGQUIT, signal_handler)
    signal.signal(signal.SIGSEGV, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)
    signal.signal(signal.SIGINT, signal_handler)
    env = os.environ.copy()
    # VESTA does not require the XLIB_SKIP_ARGB_VISUALS workaround that
    # AtomEye needed; this block is intentionally omitted.
    # Try to run the viewer:
    try:
        convert_structure_file(pd)
        spawnargs = (pd["viewer"], pd["viewer"], pd["tmpfile"], env)
        if pd["watch"]:
            signal.signal(signal.SIGCHLD, signal_handler)
            os.spawnlpe(os.P_NOWAIT, *spawnargs)
            watch_structure_file(pd)
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

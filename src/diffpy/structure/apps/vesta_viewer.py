#!/usr/bin/env python
##############################################################################
#
# diffpy.structure  by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2026 University of California, Santa Barbara.
#                   All rights reserved.
#
# File coded by:    Simon J. L. Billinge, Rundong Hua
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

import os
import re
import signal
import sys
from pathlib import Path

from diffpy.structure.structureerrors import StructureFormatError

pd = {
    "formula": None,
    "watch": False,
    "viewer": "vesta",
    "formats": ["vesta", "cif"],
}


def usage(style=None):
    """Show usage info. for ``style=="brief"`` show only first 2 lines.

    Parameters
    ----------
    style : str, optional
        The usage display style.
    """
    myname = Path(sys.argv[0]).name
    msg = __doc__.replace("vestaview", myname)
    if style == "brief":
        msg = f"{msg.splitlines()[1]}\n" f"Try `{myname} --help' for more information."
    else:
        from diffpy.structure.parsers import input_formats

        fmts = [fmt for fmt in input_formats() if fmt != "auto"]
        msg = msg.replace("inputFormats", " ".join(fmts))
    print(msg)


def version():
    """Print the script version."""
    from diffpy.structure import __version__

    print(f"vestaview {__version__}")


def load_structure_file(filename, format="auto"):
    """Load structure from the specified file.

    Parameters
    ----------
    filename : str or Path
        The path to the structure file.
    format : str, optional
        The file format, by default ``"auto"``.

    Returns
    -------
    tuple
        The loaded ``(Structure, fileformat)`` pair.
    """
    from diffpy.structure import Structure

    stru = Structure()
    parser = stru.read(str(filename), format)
    return stru, parser.format


def convert_structure_file(pd):
    """Convert ``strufile`` to a temporary file understood by the
    viewer.

    On the first call, a temporary directory is created and stored in
    ``pd``. Subsequent calls in watch mode reuse the directory.

    The VESTA viewer natively reads ``.vesta`` and ``.cif`` files, so if
    the source is already in one of the formats listed in
    ``pd["formats"]`` and no formula override is requested, the file is
    copied unchanged. Otherwise the structure is loaded and re-written in
    the first format listed in ``pd["formats"]``.

    Parameters
    ----------
    pd : dict
        The parameter dictionary containing at minimum ``"strufile"``
        and ``"formats"`` keys. It is modified in place to add
        ``"tmpdir"`` and ``"tmpfile"`` on the first call.
    """
    if "tmpdir" not in pd:
        from tempfile import mkdtemp

        pd["tmpdir"] = Path(mkdtemp())
    strufile = Path(pd["strufile"])
    tmpfile = pd["tmpdir"] / strufile.name
    tmpfile_tmp = Path(f"{tmpfile}.tmp")
    pd["tmpfile"] = tmpfile
    stru = None
    fmt = pd.get("fmt", "auto")
    if fmt == "auto":
        stru, fmt = load_structure_file(strufile)
        pd["fmt"] = fmt
    if fmt in pd["formats"] and pd["formula"] is None:
        import shutil

        shutil.copyfile(strufile, tmpfile_tmp)
        tmpfile_tmp.replace(tmpfile)
        return
    if stru is None:
        stru = load_structure_file(strufile, fmt)[0]
    if pd["formula"]:
        formula = pd["formula"]
        if len(formula) != len(stru):
            emsg = f"Formula has {len(formula)} atoms while structure has " f"{len(stru)}"
            raise RuntimeError(emsg)
        for atom, element in zip(stru, formula):
            atom.element = element
    elif fmt == "rawxyz":
        for atom in stru:
            if atom.element == "":
                atom.element = "C"
    stru.write(str(tmpfile_tmp), pd["formats"][0])
    tmpfile_tmp.replace(tmpfile)


def watch_structure_file(pd):
    """Watch ``strufile`` for modifications and reconvert when changed.

    Polls the modification timestamps of ``pd["strufile"]`` and
    ``pd["tmpfile"]`` once per second. When the source is newer, the
    file is reconverted via :func:`convert_structure_file`.

    Parameters
    ----------
    pd : dict
        The parameter dictionary as used by
        :func:`convert_structure_file`.
    """
    from time import sleep

    strufile = Path(pd["strufile"])
    tmpfile = Path(pd["tmpfile"])
    while pd["watch"]:
        if tmpfile.stat().st_mtime < strufile.stat().st_mtime:
            convert_structure_file(pd)
        sleep(1)


def clean_up(pd):
    """Remove temporary file and directory created during conversion.

    Parameters
    ----------
    pd : dict
        The parameter dictionary that may contain ``"tmpfile"`` and
        ``"tmpdir"`` entries to be removed.
    """
    tmpfile = pd.pop("tmpfile", None)
    if tmpfile is not None and Path(tmpfile).exists():
        Path(tmpfile).unlink()
    tmpdir = pd.pop("tmpdir", None)
    if tmpdir is not None and Path(tmpdir).exists():
        Path(tmpdir).rmdir()


def parse_formula(formula):
    """Parse chemical formula and return a list of elements.

    Parameters
    ----------
    formula : str
        The chemical formula string such as ``"Na4Cl4"`` or ``"H2O"``.

    Returns
    -------
    list of str
        The ordered list of element symbols with repetition matching the
        formula.

    Raises
    ------
    RuntimeError
        Raised when ``formula`` does not start with an uppercase letter
        or contains a non-integer count.
    """
    formula = re.sub(r"\s", "", formula)
    if not re.match(r"^[A-Z]", formula):
        raise RuntimeError(f"InvalidFormula '{formula}'")

    elcnt = re.split(r"([A-Z][a-z]?)", formula)[1:]
    ellst = []
    try:
        for i in range(0, len(elcnt), 2):
            element = elcnt[i]
            count = int(elcnt[i + 1]) if elcnt[i + 1] else 1
            ellst.extend([element] * count)
    except ValueError:
        emsg = f"Invalid formula, {elcnt[i + 1]!r} is not valid count"
        raise RuntimeError(emsg)
    return ellst


def die(exit_status=0, pd=None):
    """Clean up temporary files and exit with ``exit_status``.

    Parameters
    ----------
    exit_status : int, optional
        The exit code passed to :func:`sys.exit`, by default 0.
    pd : dict, optional
        The parameter dictionary forwarded to :func:`clean_up`.
    """
    clean_up({} if pd is None else pd)
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
        The signal number.
    stackframe : frame
        The current stack frame. Unused.
    """
    del stackframe
    signal.signal(signum, signal.SIG_DFL)
    if signum == signal.SIGCHLD:
        _, exit_status = os.wait()
        exit_status = (exit_status >> 8) + (exit_status & 0x00FF)
        die(exit_status, pd)
    else:
        die(1, pd)


def main():
    """Entry point for the ``vestaview`` command-line tool."""
    import getopt

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

    for option, argument in opts:
        if option in ("-f", "--formula"):
            try:
                pd["formula"] = parse_formula(argument)
            except RuntimeError as err:
                print(err, file=sys.stderr)
                die(2)
        elif option in ("-w", "--watch"):
            pd["watch"] = True
        elif option == "--viewer":
            pd["viewer"] = argument
        elif option == "--formats":
            pd["formats"] = [word.strip() for word in argument.split(",")]
        elif option in ("-h", "--help"):
            usage()
            die()
        elif option in ("-V", "--version"):
            version()
            die()
    if len(args) < 1:
        usage("brief")
        die()
    if len(args) > 1:
        print("too many structure files", file=sys.stderr)
        die(2)
    pd["strufile"] = Path(args[0])
    signal.signal(signal.SIGHUP, signal_handler)
    signal.signal(signal.SIGQUIT, signal_handler)
    signal.signal(signal.SIGSEGV, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)
    signal.signal(signal.SIGINT, signal_handler)
    env = os.environ.copy()
    try:
        convert_structure_file(pd)
        spawnargs = (
            pd["viewer"],
            pd["viewer"],
            str(pd["tmpfile"]),
            env,
        )
        if pd["watch"]:
            signal.signal(signal.SIGCHLD, signal_handler)
            os.spawnlpe(os.P_NOWAIT, *spawnargs)
            watch_structure_file(pd)
        else:
            status = os.spawnlpe(os.P_WAIT, *spawnargs)
            die(status, pd)
    except IOError as err:
        print(f"{args[0]}: {err.strerror}", file=sys.stderr)
        die(1, pd)
    except StructureFormatError as err:
        print(f"{args[0]}: {err}", file=sys.stderr)
        die(1, pd)


if __name__ == "__main__":
    main()

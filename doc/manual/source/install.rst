Installation
========================================================================

.. index:: Requirements

Requirements
------------------------------------------------------------------------

The diffpy.Structure requires Python 2.5 and the following external software:

    setuptools  -- software distribution tools for Python
    numpy       -- numerical mathematics and fast array operations for Python

On Ubuntu Linux the required software can be easily installed using
the system package manager:

    sudo aptitude install \
        python-setuptools python-numpy

For other Linux distributions use their respective package manager; note
the packages may have slightly different names.  diffpy.Structure should work
on other Unix-like operating systems and on Mac as well.  Please, search the
web for instructions how to install external dependencies on your particular
system.

.. index:: Installation

Installation
------------------------------------------------------------------------

To install the diffpy.Structure package:

    python setup.py install

By default the files are installed in the system directories, which are
usually only writeable by the root.  See the usage info
"./setup.py install --help" for options to install as a normal user under
different location.  Note that installation to non-standard directories may
require adjustments to the PATH and PYTHONPATH environment variables.

The Python setuptools library provides "easy_install" script, which can
be used to update diffpy.Structure installation or even to perform a new
install without explicit need to download and unzip the code:

    easy_install -U diffpy.Structure

This checks the package repository at http://www.diffpy.org/packages/
for any newer releases of diffpy.Structure and if present, it updates the
installation.

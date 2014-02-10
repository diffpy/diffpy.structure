diffpy.Structure - simple storage and manipulation of crystal structures

The diffpy.Structure package provides simple crystal structure operation
such as calculations of distances or angles from fractional coordinates.
It also provides imports and exports in several structure file formats
and space group definitions and utilities for symmetry expansion and 
generation of symmetry constraints on positions and thermal factors.

To learn more about diffpy.Structure library, see the examples directory
included in this distribution or the API documentation at

    http://diffpy.github.io/doc/Structure/


REQUIREMENTS

The diffpy.Structure requires Python 2.6 or 2.7 and the following software:

    setuptools  -- software distribution tools for Python
    numpy       -- numerical mathematics and fast array operations for Python

On Ubuntu Linux the required software can be easily installed using
the system package manager:

    sudo aptitude install \
        python-setuptools python-numpy
        
For Mac OS X machine with the MacPorts package manager one could do

    sudo port install \
        python27 py27-setuptools py27-numpy

When installing with MacPorts, make sure the MacPorts bin directory is the
first in the system PATH and that python27 is selected as the default
Python version in MacPorts:

    sudo port select --set python python27
    
For other Linux distributions use their respective package manager; note
the packages may have slightly different names. diffpy.Structure should work
on other Unix-like operating systems as well.  Please, search the
web for instructions how to install external dependencies on your particular
system.


INSTALLATION

The easiest option is to use the latest DiffPy-CMI release bundle from
http://www.diffpy.org/, which comes with diffpy.Structure and all other
dependencies included.

To install the diffpy.Structure package:

    python setup.py install

By default the files are installed in the system directories, which are
usually only writeable by the root.  See the usage info 
"./setup.py install --help" for options to install as a normal user under
different location.  Note that installation to non-standard directories may
require adjustments to the PATH and PYTHONPATH environment variables.

The installation integrity can be verified by changing to
the HOME directory and running

    python -m diffpy.Structure.tests.run

DEVELOPMENT and CONTRIBUTION

diffpy.Structure is an open-source software developed as a part of the
DiffPy-CMI complex modeling initiative at the Brookhaven National
Laboratory.  The diffpy.Structure sources are hosted at

    https://github.com/diffpy/diffpy.Structure

Feel free to fork the project and contribute.  To install diffpy.Structure
in a development mode, where the sources are directly used by Python
rather than copied to a system directory, use

    python setup.py develop --user


CONTACTS

For more information on diffpy.Structure please visit the project web-page:

    http://www.diffpy.org/

or email Prof. Simon Billinge at sb2896@columbia.edu

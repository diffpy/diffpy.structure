diffpy.Structure
========================================================================

storage and manipulation of crystal structure data

The diffpy.Structure package provides objects for storing atomic
coordinates, displacement parameters and other crystal structure data.
diffpy.Structure supports import and export of structure data in several
structure formats such as CIF, PDB, xyz.  It provides conversion
between fractional and absolute Cartesian coordinates, functions for
symmetry expansion from asymmetric unit and generation of symmetry
constraints for atom positions and displacement parameters.  diffpy.Structure
includes definitions of all space groups in over 500 symmetry settings.

To learn more about diffpy.Structure library, see the
`examples <https://github.com/diffpy/diffpy.Structure/tree/master/examples>`__
in the source repository and the user manual at
http://diffpy.github.io/diffpy.Structure.


REQUIREMENTS
------------------------------------------------------------------------

The diffpy.Structure requires Python 2.6 or 2.7 and the following software:

* ``setuptools`` - software distribution tools for Python
* ``NumPy`` - numerical mathematics and fast array operations for Python

On Ubuntu Linux the required software can be easily installed using
the system package manager::

   sudo aptitude install python-setuptools python-numpy

For Mac OS X machine with the MacPorts package manager the installation
command is ::

   sudo port install python27 py27-setuptools py27-numpy

When installing with MacPorts, make sure the MacPorts bin directory is the
first in the system PATH and that python27 is selected as the default
Python version in MacPorts::

   sudo port select --set python python27

For other Linux distributions use their respective package manager; note
the packages may have slightly different names.  diffpy.Structure depends
also on the `PyCifRW <http://pycifrw.berlios.de>`_ library, which should
get automatically deployed as a part of the installation process.


INSTALLATION
------------------------------------------------------------------------

Use ``easy_install`` to download and install the latest release from
`Python Package Index <https://pypi.python.org>`_ ::

   sudo easy_install diffpy.Structure

diffpy.Structure is also included in the DiffPy-CMI release bundle of
Python libraries for structure analysis at http://www.diffpy.org.

If you prefer to install from sources, make sure all required software
packages are in place and then run ::

   sudo python setup.py install

By default the files are installed in the system directories, which are
only writeable by the root user.  See the usage info from
``./setup.py install --help`` for options to install as a normal user under
different location.  Note that installation to non-standard directories may
require adjustments to the PATH and PYTHONPATH environment variables.
The installation integrity can be verified by changing to the HOME
directory and running ::

   python -m diffpy.Structure.tests.run


DEVELOPMENT
------------------------------------------------------------------------

diffpy.Structure is an open-source software developed as a part of the
DiffPy-CMI complex modeling initiative at the Brookhaven National
Laboratory.  The diffpy.Structure sources are hosted at
https://github.com/diffpy/diffpy.Structure.

Feel free to fork the project and contribute.  To install diffpy.Structure
in a development mode, where the sources are directly used by Python
rather than copied to a system directory, use ::

   python setup.py develop --user


ACKNOWLEDGEMENT
------------------------------------------------------------------------

Space group codes in *spacegroupmod.py* and *mmlibspacegroups.py*
originate from the pymmlib project, http://pymmlib.sourceforge.net.


CONTACTS
------------------------------------------------------------------------

For more information on diffpy.Structure please visit the project web-page

http://www.diffpy.org/

or email Prof. Simon Billinge at sb2896@columbia.edu.

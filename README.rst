.. image:: https://travis-ci.org/diffpy/diffpy.structure.svg?branch=master
   :target: https://travis-ci.org/diffpy/diffpy.structure

.. image:: https://codecov.io/gh/diffpy/diffpy.structure/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/diffpy/diffpy.structure


diffpy.Structure
========================================================================

storage and manipulation of crystal structure data

The diffpy.Structure package provides objects for storing atomic
coordinates, displacement parameters and other crystal structure data.
diffpy.Structure supports import and export of structure data in several
structure formats such as CIF, PDB, and xyz.  It provides conversion
between fractional and absolute Cartesian coordinates, functions for
symmetry expansion of atom sites in the asymmetric unit and generation
of symmetry constraints for atom positions and displacement parameters.
diffpy.Structure includes definitions of all space groups in over 500
symmetry settings.

To learn more about diffpy.Structure library see the
user manual at http://diffpy.github.io/diffpy.structure.


REQUIREMENTS
------------------------------------------------------------------------

The diffpy.Structure requires Python 2.7 and the following software:

* ``setuptools`` - software distribution tools for Python
* ``NumPy`` - numerical mathematics and fast array operations for Python

We recommend to use `Anaconda Python <https://www.continuum.io/downloads>`_
as it allows to install all software dependencies together with
diffpy.Structure.  For other Python distributions it is necessary to
install the required software separately.  As an example on Ubuntu
Linux the required software can be installed with ::

   sudo aptitude install python-setuptools python-numpy

diffpy.Structure also uses the
`PyCifRW <https://bitbucket.org/jamesrhester/pycifrw>`_
library, which is automatically deployed during the
installation process.


INSTALLATION
------------------------------------------------------------------------

The preferred method is to use Anaconda Python and install from the
"diffpy" channel of Anaconda packages ::

   conda config --add channels diffpy
   conda install diffpy.structure

diffpy.Structure is also included in the "diffpy-cmi" collection
of packages for structure analysis ::

   conda install diffpy-cmi

Another installation option is to use ``easy_install`` to download and
install the latest release from
`Python Package Index <https://pypi.python.org>`_ ::

   easy_install diffpy.Structure

If you prefer to install from sources, navigate to the source archive
directory and run ::

   python setup.py install

You may need to use ``sudo`` with system Python so it is allowed to
copy files to system directories.  If sudo is not available, check
the usage info from ``python setup.py install --help`` for options to
install to user-writable locations.  The installation integrity can be
verified by changing to the HOME directory and running ::

   python -m diffpy.Structure.tests.run


DEVELOPMENT
------------------------------------------------------------------------

diffpy.Structure is an open-source software developed as a part of the
DiffPy-CMI complex modeling initiative at the Brookhaven National
Laboratory.  The diffpy.Structure sources are hosted at
https://github.com/diffpy/diffpy.structure.

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

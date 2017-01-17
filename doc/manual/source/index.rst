.. diffpy.structure documentation master file, created by
   sphinx-quickstart on Tue Oct 22 12:02:48 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

####################################################
diffpy.structure documentation
####################################################

Software version |release|.

Last updated |today|.

diffpy.structure - simple storage and manipulation of crystal structures

The diffpy.structure package provides objects for storing atomic
coordinates, displacement parameters and other crystal structure data.
diffpy.structure supports import and export of structure data in several
structure formats such as CIF, PDB, xyz.  It provides conversion
between fractional and absolute Cartesian coordinates, functions for
symmetry expansion from asymmetric unit and generation of symmetry
constraints for atom positions and displacement parameters.  diffpy.structure
includes definitions of all space groups in over 500 symmetry settings.

To learn more about diffpy.structure library, see the examples directory
included in this distribution or the API documentation

===================
Disclaimer
===================

.. literalinclude:: ../../../LICENSE.txt

================
Acknowledgments
================

Developers
-----------

diffpy.structure is developed and maintained by

.. literalinclude:: ../../../AUTHORS.txt


Space group codes in *spacegroupmod.py* and *mmlibspacegroups.py*
originate from the pymmlib project, http://pymmlib.sourceforge.net.


======================================
Installation
======================================

See the `README.rst <https://github.com/diffpy/diffpy.structure#requirements>`_
file included with the distribution.

======================================
API and Index
======================================

.. toctree::
   :maxdepth: 2

   api/diffpy.structure.rst

* :ref:`genindex`

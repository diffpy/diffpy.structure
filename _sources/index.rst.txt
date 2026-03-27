#######
|title|
#######

.. |title| replace:: diffpy.structure documentation

``diffpy.structure`` - Crystal structure container and parsers for structure formats.

| Software version |release|
| Last updated |today|.

===============
Getting started
===============

The diffpy.structure package provides objects for storing atomic
coordinates, displacement parameters and other crystal structure data.
diffpy.structure supports import and export of structure data in several
structure formats such as CIF, PDB, xyz.  It provides conversion
between fractional and absolute Cartesian coordinates, functions for
symmetry expansion from asymmetric unit and generation of symmetry
constraints for atom positions and displacement parameters.  diffpy.structure
includes definitions of all space groups in over 500 symmetry settings.

=======
Authors
=======

``diffpy.structure`` is developed by Chris Farrow, Pavol Juhas, Simon Billinge, Billinge Group members. This project is maintained by Simon Billinge. For a detailed list of contributors see
https://github.com/diffpy/diffpy.structure/graphs/contributors.

============
Installation
============

See the `README <https://github.com/diffpy/diffpy.structure#installation>`_
file included with the distribution.

================
Acknowledgements
================

Space group codes in *spacegroupmod.py* and *mmlibspacegroups.py*
originate from the pymmlib project, http://pymmlib.sourceforge.net.
Less common settings of space groups were generating using the
Computational Crystallography Toolbox,
http://cctbx.sourceforge.net.

``diffpy.structure`` is built and maintained with `scikit-package <https://scikit-package.github.io/scikit-package/>`_.

.. index:: citation, reference

=========
Reference
=========

If you use this program for a scientific research that leads
to publication, we ask that you acknowledge use of the program
by citing the following paper in your publication:

   P. Juhás, C. L. Farrow, X. Yang, K. R. Knox and S. J. L. Billinge,
   `Complex modeling: a strategy and software program for combining
   multiple information sources to solve ill posed structure and
   nanostructure inverse problems
   <http://dx.doi.org/10.1107/S2053273315014473>`__,
   *Acta Crystallogr. A* **71**, 562-568 (2015).

=================
Table of contents
=================
.. toctree::
   :maxdepth: 2

   Package API <api/diffpy.structure>
   release
   license

=======
Indices
=======

* :ref:`genindex`
* :ref:`search`

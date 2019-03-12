:tocdepth: 2

Package API
###########


Submodules
==========

.. toctree::
    :titlesonly:

    mod_lattice
    diffpy.structure.parsers
    diffpy.structure.applications
    diffpy.structure.expansion


diffpy.structure
================

.. module:: diffpy.structure

The top-level diffpy.structure module provides the following objects.


Functions
---------

.. autofunction:: loadStructure


Classes
-------

:py:class:`~atom.Atom`
   describes one atom site in the structure - its type, position,
   displacement parameters and so forth.

:py:class:`~lattice.Lattice`
   depicts general coordinate system and associated operations.

:py:class:`~structure.Structure`
   collection of :class:`!Atom` objects that form the structure together
   with associated :py:class:`!Lattice`.


Exceptions
----------

.. autoclass:: StructureFormatError

.. autoclass:: LatticeError

.. autoclass:: SymmetryError


diffpy.structure.spacegroups
============================

.. automodule:: diffpy.structure.spacegroups
    :members:
    :undoc-members:
    :show-inheritance:

diffpy.structure.structureerrors
================================

.. automodule:: diffpy.structure.structureerrors
    :members:
    :undoc-members:
    :show-inheritance:

diffpy.structure.symmetryutilities
==================================

.. automodule:: diffpy.structure.symmetryutilities
    :members:
    :undoc-members:
    :show-inheritance:

diffpy.structure.atom
=====================

.. automodule:: diffpy.structure.atom
    :members:
    :show-inheritance:

diffpy.structure.mmlibspacegroups
=================================

.. automodule:: diffpy.structure.mmlibspacegroups
    :members:
    :undoc-members:
    :show-inheritance:

diffpy.structure.pdffitstructure
================================

.. automodule:: diffpy.structure.pdffitstructure
    :members:
    :undoc-members:
    :show-inheritance:

diffpy.structure.sgtbxspacegroups
=================================

.. automodule:: diffpy.structure.sgtbxspacegroups
    :members:
    :undoc-members:
    :show-inheritance:

diffpy.structure.spacegroupmod
==============================

.. automodule:: diffpy.structure.spacegroupmod
    :members:
    :undoc-members:
    :special-members: __call__, __eq__

diffpy.structure.structure
==========================

.. automodule:: diffpy.structure.structure
    :members:
    :undoc-members:
    :show-inheritance:

diffpy.structure.utils
======================

.. automodule:: diffpy.structure.utils
    :members:
    :undoc-members:
    :show-inheritance:

diffpy.structure.version
========================

.. automodule:: diffpy.structure.version
    :members:
    :undoc-members:
    :show-inheritance:

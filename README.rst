|Icon| |title|_
===============

.. |title| replace:: diffpy.structure
.. _title: https://diffpy.github.io/diffpy.structure

.. |Icon| image:: https://avatars.githubusercontent.com/diffpy
        :target: https://diffpy.github.io/diffpy.structure
        :height: 100px

|PyPi| |Forge| |PythonVersion| |PR|

|CI| |Codecov| |Black| |Tracking|

.. |Black| image:: https://img.shields.io/badge/code_style-black-black
        :target: https://github.com/psf/black

.. |CI| image:: https://github.com/diffpy/diffpy.structure/actions/workflows/matrix-and-codecov-on-merge-to-main.yml/badge.svg
        :target: https://github.com/diffpy/diffpy.structure/actions/workflows/matrix-and-codecov-on-merge-to-main.yml

.. |Codecov| image:: https://codecov.io/gh/diffpy/diffpy.structure/branch/main/graph/badge.svg
        :target: https://codecov.io/gh/diffpy/diffpy.structure

.. |Forge| image:: https://img.shields.io/conda/vn/conda-forge/diffpy.structure
        :target: https://anaconda.org/conda-forge/diffpy.structure

.. |PR| image:: https://img.shields.io/badge/PR-Welcome-29ab47ff

.. |PyPi| image:: https://img.shields.io/pypi/v/diffpy.structure
        :target: https://pypi.org/project/diffpy.structure/

.. |PythonVersion| image:: https://img.shields.io/pypi/pyversions/diffpy.structure
        :target: https://pypi.org/project/diffpy.structure/

.. |Tracking| image:: https://img.shields.io/badge/issue_tracking-github-blue
        :target: https://github.com/diffpy/diffpy.structure/issues

Crystal structure container and parsers for structure formats.

The diffpy.structure package provides objects for storing atomic
coordinates, displacement parameters and other crystal structure data.
diffpy.structure supports import and export of structure data in several
structure formats such as CIF, PDB, and xyz.  It provides conversion
between fractional and absolute Cartesian coordinates, functions for
symmetry expansion of atom sites in the asymmetric unit and generation
of symmetry constraints for atom positions and displacement parameters.
diffpy.structure includes definitions of all space groups in over 500
symmetry settings.


For more information about the diffpy.structure library, please consult our `online documentation <https://diffpy.github.io/diffpy.structure>`_.

Citation
--------

If you use this program for a scientific research that leads
to publication, we ask that you acknowledge use of the program
by citing the following paper in your publication:

   P. Juhás, C. L. Farrow, X. Yang, K. R. Knox and S. J. L. Billinge,
   `Complex modeling: a strategy and software program for combining
   multiple information sources to solve ill posed structure and
   nanostructure inverse problems
   <http://dx.doi.org/10.1107/S2053273315014473>`__,
   *Acta Crystallogr. A* **71**, 562-568 (2015).

Installation
------------

The preferred method is to use `Miniconda Python
<https://docs.conda.io/projects/miniconda/en/latest/miniconda-install.html>`_
and install from the "conda-forge" channel of Conda packages.

To add "conda-forge" to the conda channels, run the following in a terminal. ::

        conda config --add channels conda-forge

We want to install our packages in a suitable conda environment.
The following creates and activates a new environment named ``diffpy.structure_env`` ::

        conda create -n diffpy.structure_env python=3
        conda activate diffpy.structure_env

Then, to fully install ``diffpy.structure`` in our active environment, run ::

        conda install diffpy.structure

Another option is to use ``pip`` to download and install the latest release from
`Python Package Index <https://pypi.python.org>`_.
To install using ``pip`` into your ``diffpy.structure_env`` environment, type ::

        pip install diffpy.structure

If you prefer to install from sources, after installing the dependencies, obtain the source archive from
`GitHub <https://github.com/diffpy/diffpy.structure/>`_. Once installed, ``cd`` into your ``diffpy.structure`` directory
and run the following ::

        pip install .

Support and Contribute
----------------------

`Diffpy user group <https://groups.google.com/g/diffpy-users>`_ is the discussion forum for general questions and discussions about the use of diffpy.structure. Please join the diffpy.structure users community by joining the Google group. The diffpy.structure project welcomes your expertise and enthusiasm!

If you see a bug or want to request a feature, please `report it as an issue <https://github.com/diffpy/diffpy.structure/issues>`_ and/or `submit a fix as a PR <https://github.com/diffpy/diffpy.structure/pulls>`_. You can also post it to the `Diffpy user group <https://groups.google.com/g/diffpy-users>`_.

Feel free to fork the project and contribute. To install diffpy.structure
in a development mode, with its sources being directly used by Python
rather than copied to a package directory, use the following in the root
directory ::

        pip install -e .

To ensure code quality and to prevent accidental commits into the default branch, please set up the use of our pre-commit
hooks.

1. Install pre-commit in your working environment by running ``conda install pre-commit``.

2. Initialize pre-commit (one time only) ``pre-commit install``.

Thereafter your code will be linted by black and isort and checked against flake8 before you can commit.
If it fails by black or isort, just rerun and it should pass (black and isort will modify the files so should
pass after they are modified). If the flake8 test fails please see the error messages and fix them manually before
trying to commit again.

Improvements and fixes are always appreciated.

Before contribuing, please read our `Code of Conduct <https://github.com/diffpy/diffpy.structure/blob/main/CODE_OF_CONDUCT.rst>`_.

Acknowledgement
---------------

Space group codes in *spacegroupmod.py* and *mmlibspacegroups.py*
originate from the `pymmlib project <http://pymmlib.sourceforge.net>`_.

Less common settings of space groups were generating using the
`Computational Crystallography Toolbox <http://cctbx.sourceforge.net>`_.

Contact
-------

For more information on diffpy.structure please visit the project `web-page <https://diffpy.github.io/>`_ or email Prof. Simon Billinge at sb2896@columbia.edu.

=============
Release notes
=============

.. current developments

3.4.0
=====

**Added:**

* Added ``set_latt_parms`` method into ``Lattice`` class
* Added ``set_new_latt_base_vec`` method into ``Lattice`` class
* Added ``parse_lines`` method in ``p_auto.py``
* Added ``parse_lines`` method in ``p_cif.py``
* Added ``parse_lines`` method in ``p_discus.py``
* Added ``parse_lines`` method in ``p_pdb.py``
* Added ``parse_lines`` method in ``p_pdffit.py``
* Added ``parse_lines`` method in ``p_rawxyz.py``
* Added ``parse_lines`` method in ``p_xcfg.py``
* Added ``parse_lines`` method in ``p_xyz.py``
* Added ``parse_lines`` method in ``structureparser.py``
* Added ``_suppress_cif_parser_output`` method in ``p_cif.py``
* Add deprecation warning for ``diffpy.Structure`` import.
* Added `diffpy.structure.Structure.add_new_atom` in replace of `addNewAtom`
* Added ``load_structure_file`` method in ``apps/anyeye.py``
* Added ``convert_structure_file`` method in ``apps/anyeye.py``
* Added ``watch_structure_file`` method in ``apps/anyeye.py``
* Added ``clean_up`` method in ``apps/anyeye.py``
* Added ``parse_formula`` method in ``apps/anyeye.py``
* Added ``signal_handler`` method in ``apps/anyeye.py``
* Added method ``load_structure`` in ``__init__.py``
* Added `diffpy.structure.Structure.assign_unique_labels` in replace of `assignUniqueLabels`
* Support for Python 3.14
* Added ``place_in_lattice`` method to ``Structure``
* Added ``read_structure`` method to ``Structure``
* Added ``write_structure`` method to ``Structure``
* Added ``position_formula`` method in ``GeneratorSite`` class
* Added ``u_formula`` method in ``GeneratorSite`` class
* Added ``eq_index`` method in ``GeneratorSite`` class
* Added ``prune_formula_dictionary`` method in ``symmetryutilities.py``
* Added ``_link_atom_attribute`` method in ``diffpy.structure.utils``
* Added ``msd_latt`` method in ``atom.py``
* Added ``msd_cart`` method in ``atom.py``
* Added ``_get_uij`` method in ``atom.py``
* Added ``_set_uij`` method in ``atom.py``
* Added ``parse_file`` method in ``structureparser.py``
* Added ``parse_lines`` method in ``p_cif.py``
* Added ``parse_lines`` method in ``p_auto.py``
* Added parser for vesta specific files and viewer for vesta
* Added ``atom_bare_symbol`` method in ``utils.py``
* Added ``_get_ordered_formats`` method in ``p_auto.py``
* Added ``_wrap_parse_method`` method in ``p_auto.py``
* Added ``_tr_atom_site_u_iso_or_equiv`` method in ``p_cif.py``
* Added ``_tr_atom_site_b_iso_or_equiv`` method in ``p_cif.py``
* Added ``_tr_atom_site_aniso_u_11`` method in ``p_cif.py``
* Added ``_tr_atom_site_aniso_u_22`` method in ``p_cif.py``
* Added ``_tr_atom_site_aniso_u_33`` method in ``p_cif.py``
* Added ``_tr_atom_site_aniso_u_12`` method in ``p_cif.py``
* Added ``_tr_atom_site_aniso_u_13`` method in ``p_cif.py``
* Added ``_tr_atom_site_aniso_u_23`` method in ``p_cif.py``
* Added ``_tr_atom_site_aniso_b_11`` method in ``p_cif.py``
* Added ``_tr_atom_site_aniso_b_22`` method in ``p_cif.py``
* Added ``_tr_atom_site_aniso_b_33`` method in ``p_cif.py``
* Added ``_tr_atom_site_aniso_b_12`` method in ``p_cif.py``
* Added ``_tr_atom_site_aniso_b_13`` method in ``p_cif.py``
* Added ``_tr_atom_site_aniso_b_23`` method in ``p_cif.py``
* Added ``get_symop`` method in ``parsers/p_cif.py``
* Added ``get_space_group`` method in ``spacegroups.py``
* Added ``find_space_group`` method in ``spacegroups.py``
* Added ``is_space_group_identifier`` method in ``spacegroups.py``
* Added ``_hash_symop_list`` method in ``spacegroups.py``
* Added ``_build_sg_lookup_table`` method in ``spacegroups.py``
* Added ``_get_sg_hash_lookup_table`` method in ``spacegroups.py``
* Added ``read_structure`` method into ``PDFFitStructure`` class
* Added ``cell_parms`` method into ``Lattice`` class
* Added ``_find_constraints`` method in ``SymmetryConstraints`` class
* Added ``pos_parm_symbols`` method in ``SymmetryConstraints`` class
* Added ``pos_parm_values`` method in ``SymmetryConstraints`` class
* Added ``u_parm_symbols`` method in ``SymmetryConstraints`` class
* Added ``u_parm_values`` method in ``SymmetryConstraints`` class
* Added ``u_formulas`` method in ``SymmetryConstraints`` class
* Added `diffpy.structure.Structure.get_last_atom` in replace of `getLastAtom`
* Added ``get_parser`` method in ``p_auto.py``
* Added ``get_parser`` method in ``p_cif.py``
* Added ``get_parser`` method in ``p_discus.py``
* Added ``get_parser`` method in ``p_pdb.py``
* Added ``get_parser`` method in ``p_pdffit.py``
* Added ``get_parser`` method in ``p_rawxyz.py``
* Added ``get_parser`` method in ``p_xcfg.py``
* Added ``get_parser`` method in ``p_xyz.py``
* Added ``get_parser`` method in ``parsers/__init__.py``
* Added ``position_formulas`` method in ``SymmetryConstraints`` class
* Added ``position_formulas_pruned`` method in ``SymmetryConstraints`` class
* Added ``u_formulas_pruned`` method in ``SymmetryConstraints`` class
* Added ``_parse_cif_data_source`` method in ``p_cif.py``
* Added ``_parse_cif_block`` method in ``p_cif.py``
* Added ``to_lines`` method in ``p_cif.py``
* Added ``to_lines`` method in ``p_pdb.py``
* Added ``to_lines`` method in ``p_rawxyz.py``
* Added ``to_lines`` method in ``p_xcfg.py``
* Added ``to_lines`` method in ``p_xyz.py``
* Added ``to_lines`` method in ``structureparser.py``
* Added ``_lines_iterator`` method in ``p_discus.py``
* Added ``to_lines`` method in ``p_discus.py``
* Added ``is_space_group_latt_parms`` method in ``symmetryutilities.py``
* Added ``is_constant_formula`` method in ``symmetryutilities.py``
* Added ``find_center`` method in ``expansion/shapeutils.py``
* Added ``make_sphere`` method in ``expansion/makeellipsoid.py``
* Added ``make_ellipsoid`` method in ``expansion/makeellipsoid.py``
* Added ``position_difference`` method in ``symmetryutilities.py``
* Added ``nearest_site_index`` method in ``symmetryutilities.py``
* Added ``_find_invariants`` method in ``symmetryutilities.py``
* Added ``equal_positions`` method in ``symmetryutilities.py``
* Added ``expand_position`` method in ``symmetryutilities.py``
* Added ``null_space`` method in ``symmetryutilities.py``
* Added ``input_formats`` method in ``parsers/__init__.py``
* Added ``output_formats`` method in ``parsers/__init__.py``
* Added ``title_lines`` method in ``p_pdb.py``
* Added ``cryst1_lines`` method in ``p_pdb.py``
* Added ``atom_lines`` method in ``p_pdb.py``
* Added ``convert_fp_num_to_signed_rational`` method in ``GeneratorSite`` class
* Added ``_find_null_space`` method in ``GeneratorSite`` class
* Added ``_find_pos_parameters`` method in ``GeneratorSite`` class
* Added ``_find_u_space`` method in ``GeneratorSite`` class
* Added ``_find_u_parameters`` method in ``GeneratorSite`` class
* Added ``_find_eq_uij`` method in ``GeneratorSite`` class

**Changed:**

* Changed private method ``__emptySharedStructure`` to ``__empty_shared_structure``

**Deprecated:**

* Deprecated ``setLatPar`` method in ``Lattice`` class for removal in version 4.0.0
* Deprecated ``setLatBase`` method in ``Lattice`` class for removal in version 4.0.0
* Deprecated ``parseLines`` method in ``p_auto.py`` for removal in version 4.0.0
* Deprecated ``parseLines`` method in ``p_cif.py`` for removal in version 4.0.0
* Deprecated ``parseLines`` method in ``p_discus.py`` for removal in version 4.0.0
* Deprecated ``parseLines`` method in ``p_pdb.py`` for removal in version 4.0.0
* Deprecated ``parseLines`` method in ``p_pdffit.py`` for removal in version 4.0.0
* Deprecated ``parseLines`` method in ``p_rawxyz.py`` for removal in version 4.0.0
* Deprecated ``parseLines`` method in ``p_xcfg.py`` for removal in version 4.0.0
* Deprecated ``parseLines`` method in ``p_xyz.py`` for removal in version 4.0.0
* Deprecated ``parseLines`` method in ``structureparser.py`` for removal in version 4.0.0
* Deprecated `diffpy.structure.Structure.addNewAtom` method for removal in version 4.0.0
* Deprecated ``loadStructureFile`` method in ``apps/anyeye.py`` for removal in version 4.0.0
* Deprecated ``convertStructureFile`` method in ``apps/anyeye.py`` for removal in version 4.0.0
* Deprecated ``watchStructureFile`` method in ``apps/anyeye.py`` for removal in version 4.0.0
* Deprecated ``cleanUp`` method in ``apps/anyeye.py`` for removal in version 4.0.0
* Deprecated ``parseFormula`` method in ``apps/anyeye.py`` for removal in version 4.0.0
* Deprecated ``signalHandler`` method in ``apps/anyeye.py`` for removal in version 4.0.0
* Deprecated method ``loadStructure`` in ``__init__.py`` for removal in version 4.0.0
* Deprecated `diffpy.structure.Structure.assignUniqueLabels` for removal in 4.0.0
* Deprecated ``placeInLattice`` method of ``Structure`` for removal in version 4.0.0
* Deprecated ``readStr`` method of ``Structure`` for removal in version 4.0.0
* Deprecated ``writeStr`` method of ``Structure`` for removal in version 4.0.0
* Deprecated ``positionFormula`` method in ``GeneratorSite`` class for removal in version 4.0.0
* Deprecated ``UFormula`` method in ``GeneratorSite`` class for removal in version 4.0.0
* Deprecated ``eqIndex`` method in ``GeneratorSite`` class for removal in version 4.0.0
* Deprecated ``pruneFormulaDictionary`` method in ``symmetryutilities.py`` for removal in version 4.0.0
* Derecated ``_linkAtomAttribute`` method in ``diffpy.structure.utils`` for removal in version 4.0.0
* Deprecated ``msdLat`` method in ``atom.py`` for removal in version 4.0.0
* Deprecated ``msdCart`` method in ``atom.py`` for removal in version 4.0.0
* Deprecated ``parse_file`` method in ``structureparser.py`` for removal in version 4.0.0
* Deprecated ``parse_file`` method in ``p_cif.py`` for removal in version 4.0.0
* Deprecated ``parse_file`` method in ``p_auto.py`` for removal in version 4.0.0
* Deprecated ``atomBareSymbol`` method in ``utils.py`` for removal in version 4.0.0
* Deprecated ``getSymOp`` method in ``parsers/p_cif.py`` for removal in version 4.0.0
* Deprecated ``GetSpaceGroup`` method in ``spacegroups.py`` for removal in version 4.0.0
* Deprecated ``IsSpaceGroupIdentifier`` method in ``spacegroups.py`` for removal in version 4.0.0
* Deprecated ``FindSpaceGroup`` method in ``spacegroups.py`` for removal in version 4.0.0
* Deprecated ``_hashSymOpList`` method in ``spacegroups.py`` for removal in version 4.0.0
* Deprecated ``readStr`` method in ``PDFFitStructure`` class for removal in version 4.0.0
* Deprecated ``abcABG`` method in ``Lattice`` class for removal in version 4.0.0
* Deprecated ``posparSymbols`` method in ``SymmetryConstraints`` class for removal in version 4.0.0
* Deprecated ``posparValues`` method in ``SymmetryConstraints`` class for removal in version 4.0.0
* Deprecated ``UparSymbols`` method in ``SymmetryConstraints`` class for removal in version 4.0.0
* Deprecated ``UparValues`` method in ``SymmetryConstraints`` class for removal in version 4.0.0
* Deprecated ``UFormulas`` method in ``SymmetryConstraints`` class for removal in version 4.0.0
* Deprecated `diffpy.structure.Structure.getLastAtom` for removal in version 4.0.0
* Deprecated ``getParser`` method in ``p_auto.py`` for removal in version 4.0.0
* Deprecated ``getParser`` method in ``p_cif.py`` for removal in version 4.0.0
* Deprecated ``getParser`` method in ``p_discus.py`` for removal in version 4.0.0
* Deprecated ``getParser`` method in ``p_pdb.py`` for removal in version 4.0.0
* Deprecated ``getParser`` method in ``p_pdffit.py`` for removal in version 4.0.0
* Deprecated ``getParser`` method in ``p_rawxyz.py`` for removal in version 4.0.0
* Deprecated ``getParser`` method in ``p_xcfg.py`` for removal in version 4.0.0
* Deprecated ``getParser`` method in ``p_xyz.py`` for removal in version 4.0.0
* Deprecated ``getParser`` method in ``parsers/__init__.py`` for removal in version 4.0.0
* Deprecated ``positionFormulas`` method in ``SymmetryConstraints`` class for removal in version 4.0.0
* Deprecated ``positionFormulasPruned`` method in ``SymmetryConstraints`` class for removal in version 4.0.0
* Deprecated ``UFormulasPruned`` method in ``SymmetryConstraints`` class for removal in version 4.0.0
* Deprecated ``toLines`` method in ``p_cif.py`` for removal in version 4.0.0
* Deprecated ``toLines`` method in ``p_pdb.py`` for removal in version 4.0.0
* Deprecated ``toLines`` method in ``p_rawxyz.py`` for removal in version 4.0.0
* Deprecated ``toLines`` method in ``p_xcfg.py`` for removal in version 4.0.0
* Deprecated ``toLines`` method in ``p_xyz.py`` for removal in version 4.0.0
* Deprecated ``toLines`` method in ``structureparser.py`` for removal in version 4.0.0
* Deprecated ``toLines`` method in ``p_discus.py`` for removal in version 4.0.0
* Deprecated ``isSpaceGroupLatPar`` method in ``symmetryutilities.py`` for removal in version 4.0.0
* Deprecated ``isconstantFormula`` method in ``symmetryutilities.py`` for removal in version 4.0.0
* Deprecated ``findCenter`` method in ``expansion/shapeutils.py`` for removal in version 4.0.0
* Deprecated ``makeSphere`` method in ``expansion/makeellipsoid.py`` for removal in version 4.0.0
* Deprecated ``makeEllipsoid`` method in ``expansion/makeellipsoid.py`` for removal in version 4.0.0
* Deprecated ``positionDifference`` method in ``symmetryutilities.py`` for removal in version 4.0.0
* Deprecated ``nearestSiteIndex`` method in ``symmetryutilities.py`` for removal in version 4.0.0
* Deprecated ``equalPositions`` method in ``symmetryutilities.py`` for removal in version 4.0.0
* Deprecated ``expandPosition`` method in ``symmetryutilities.py`` for removal in version 4.0.0
* Deprecated ``nullSpace`` method in ``symmetryutilities.py`` for removal in version 4.0.0
* Deprecated ``inputFormats`` method in ``parsers/__init__.py`` for removal in version 4.0.0
* Deprecated ``outputFormats`` method in ``parsers/__init__.py`` for removal in version 4.0.0
* Deprecated ``titleLines`` method in ``p_pdb.py`` for removal in version 4.0.0
* Deprecated ``crystl1Lines`` method in ``p_pdb.py`` for removal in version 4.0.0
* Deprecated ``atomLines`` method in ``p_pdb.py`` for removal in version 4.0.0
* Deprecated ``signedRatStr`` method in in ``GeneratorSite`` class for removal in version 4.0.0

**Fixed:**

* Fixed ``load_structure`` with successfully loading `Path` object
* Fix deprecation for open file for Python 3.14

**Removed:**

* Support for Python 3.11


3.3.1
=====

**Added:**

* Spelling check via Codespell in pre-commit
* Coverage report in each PR

**Changed:**

* Use the names CODE-OF-CONDUCT.rst, docs and requirements/tests.txt according to the new group standard.

**Fixed:**

* Let ``diffpy.structure`` pass the tests with ``pycifrw`` installed from ``PyPI``.


3.3.0
=====

**Added:**

* Support for Python 3.13

**Deprecated:**

* Support for Python 3.10


3.2.3
=====

**Added:**

* Use GitHub Actions to build, release, upload to PyPI
* Added issue template for PyPI/GitHub release
* Include GitHub Issues templates for bug report and feature request

**Fixed:**

* Add getting started section and re-arrange install success check instructions
* Added terminal script for transtru app in pyproject.toml
* Changed requires-python to align with classifiers


3.2.2
=====

**Added:**

* Unit test for version.py
* support for numpy >=2.0

**Fixed:**

* tests folder at the root of the repo
* Add pip dependencies under pip.txt and conda dependencies under conda.txt
* element/label itemsize to 5 in _linkAtomAttribute to support numpy >=2.0
* Recookiecut package to the group standard



3.2.1
=====



3.2.0
=====

**Changed:**

* Removed support for Python 2
* This version only supporting Python 3.10, 3.11, 3.12
* All docstrings style updated to numpydoc

**Deprecated:**

* Deprecated the `diffpy.structure.applications` module. Use
  `diffpy.structure.apps` instead

**Removed:**

* Removed all `six` compatibility code

**Fixed:**

* Repo structure modified to the new diffpy standard



Version 3.1.0 - 2022-12-04
--------------------------

**Added**

- Compatibility with Python 3.10, 3.9, 3.8

**Changed**

**Deprecated**

**Removed**

- Remove the support for Python 3.5, 3.6.

**Fixed**

Version 3.0.2 - 2022-10-12
--------------------------

**Added**

- A string representation of `SpaceGroup` with key information.

**Changed**

- Bumped minimum `PyCifRW` version to `4.4.3`.

**Deprecated**

**Removed**

**Fixed**

- Handling of paths on Windows when using the `P_cif` parser.

Version 3.0.1 - 2019-06-27
--------------------------

**Added**

- Function `FindSpaceGroup` for space group lookup from its list
  of symmetry operations.

**Changed**

- Reuse existing `SpaceGroup` instance when loading a CIF file.
- Improve check of SpaceGroup identifiers in `GetSpaceGroup`.
- When loading CIF file, preset `Atom.anisotropy` according
  to symmetry constraints at each site.  Adhere to specific
  ADP type when specified in the CIF.

**Removed**

- Unused attribute `SpaceGroup.alt_name`.

**Fixed**

- Fix inconsistent (`Atom`, `Structure`) pickle.  Preserve `Atom`
  ownership in a `Structure` after pickling and unpickling.
- Spuriously linked array-view values after `stru.xyz = 0`.
- Preserve scalar value type when setting `stru.occupancy = value`.
- Process unknown CIF occupancy "?" as an occupancy of 1.
- Incorrect `SymOp` list for spacegroup "B11m" (number 1008).


Version 3.0.0 - 2019-03-11
--------------------------

Notable differences from version 1.3.5.

**Added**

- Compatibility with Python 3.7, 3.6, 3.5 in addition to 2.7.
- Aliases for 17 non-standard space group names from cctbx.
- Support for intersphinx links to Python and NumPy documentation.
- Dependency and use of the `six` PY2/PY3 compatibility package.
- Documentation hosting at readthedocs.org.

**Changed**

- Rename the package and all its module names to lowercase.
- Use UTF-8 encoding when writing structure files.
- Refactor parsing of XCFG format.  Avoid use of generated code.
- Refactor all starred imports to explicit so they can be checked.
- Adopt napoleon style for docstrings.
- Update docstrings for `Atom`, `Lattice`, `SymOp`, `SpaceGroup`.
- Switch to platform-independent "noarch" Anaconda package.

**Deprecated**

- Old camel case module names such as `diffpy.Structure`.
- Variable `__gitsha__` in the `version` module which was renamed
  to `__git_commit__`.

**Removed**

- Unused exception `IsotropyError`.
- Unused class `BRAtomsStructure` and associated parser.

**Fixed**

- Loading of empty CIF files with no specified sites.
- Parsing of CIFs with `?` value for unknown displacement parameters.
- Symmetry constraint equations for ADPs so they avoid self-reference.
- Use `StructureFormatError` exception for CIF with unknown space group.
- Open files within the `with` context so they get closed when done.
- Invalid escape sequences in string values.

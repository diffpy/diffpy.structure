# Release notes

## Unreleased - Version 3.0.2

### Added

- A string representation of `SpaceGroup` with key information.

### Changed

### Deprecated

### Removed

### Fixed


## Version 3.0.1 – 2019-06-27

### Added

- Function `FindSpaceGroup` for space group lookup from its list
  of symmetry operations.

### Changed

- Reuse existing `SpaceGroup` instance when loading a CIF file.
- Improve check of SpaceGroup identifiers in `GetSpaceGroup`.
- When loading CIF file, preset `Atom.anisotropy` according
  to symmetry constraints at each site.  Adhere to specific
  ADP type when specified in the CIF.

### Removed

- Unused attribute `SpaceGroup.alt_name`.

### Fixed

- Fix inconsistent (`Atom`, `Structure`) pickle.  Preserve `Atom`
  ownership in a `Structure` after pickling and unpickling.
- Spuriously linked array-view values after `stru.xyz = 0`.
- Preserve scalar value type when setting `stru.occupancy = value`.
- Process unknown CIF occupancy "?" as an occupancy of 1.
- Incorrect `SymOp` list for spacegroup "B11m" (number 1008).


## Version 3.0.0 – 2019-03-11

Notable differences from version 1.3.5.

### Added

- Compatibility with Python 3.7, 3.6, 3.5 in addition to 2.7.
- Aliases for 17 non-standard space group names from cctbx.
- Support for intersphinx links to Python and NumPy documentation.
- Dependency and use of the `six` PY2/PY3 compatibility package.
- Documentation hosting at readthedocs.org.

### Changed

- Rename the package and all its module names to lowercase.
- Use UTF-8 encoding when writing structure files.
- Refactor parsing of XCFG format.  Avoid use of generated code.
- Refactor all starred imports to explicit so they can be checked.
- Adopt napoleon style for docstrings.
- Update docstrings for `Atom`, `Lattice`, `SymOp`, `SpaceGroup`.
- Switch to platform-independent "noarch" Anaconda package.

### Deprecated

- Old camel case module names such as `diffpy.Structure`.
- Variable `__gitsha__` in the `version` module which was renamed
  to `__git_commit__`.

### Removed

- Unused exception `IsotropyError`.
- Unused class `BRAtomsStructure` and associated parser.

### Fixed

- Loading of empty CIF files with no specified sites.
- Parsing of CIFs with `?` value for unknown displacement parameters.
- Symmetry constraint equations for ADPs so they avoid self-reference.
- Use `StructureFormatError` exception for CIF with unknown space group.
- Open files within the `with` context so they get closed when done.
- Invalid escape sequences in string values.

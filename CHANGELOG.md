# Release notes

## Unreleased - Version 3.0.0

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

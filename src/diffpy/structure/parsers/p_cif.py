#!/usr/bin/env python
##############################################################################
#
# diffpy.structure  by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2007 trustees of the Michigan State University.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################

"""Parser for basic CIF file format.

Attributes
----------
rx_float : re.Pattern
    Constant regular expression for `leading_float()`.
symvec : dict
    Helper dictionary for `getSymOp()`.

Note
----
References: https://www.iucr.org/resources/cif
"""

import io
import re
import sys
from contextlib import contextmanager

import numpy

from diffpy.structure import Atom, Lattice, Structure
from diffpy.structure.parsers import StructureParser
from diffpy.structure.structureerrors import StructureFormatError

# ----------------------------------------------------------------------------


class P_cif(StructureParser):
    """Simple parser for CIF structure format.

    Reads Structure from the first block containing _atom_site_label key.
    Following blocks, if any, are ignored.

    Parameters
    ----------
    eps : float, Optional
        Fractional coordinates cutoff for duplicate positions.
        When ``None`` use the default for `ExpandAsymmetricUnit`: ``1.0e-5``.

    Attributes
    ----------
    format : str
        Structure format name.
    ciffile : CifFile
        Instance of `CifFile` from `PyCifRW`.
    stru : Structure
        `Structure` instance used for CIF input or output.
    spacegroup : SpaceGroup
        Instance of `SpaceGroup` used for symmetry expansion.
    eps : float
        Resolution in fractional coordinates for non-equal positions.
        Used for expansion of asymmetric unit.
    eau : ExpandAsymmetricUnit
        Instance of `ExpandAsymmetricUnit` from `SymmetryUtilities`.
    asymmetric_unit : list
        List of `Atom` instances for the original asymmetric unit in the CIF file.
    labelindex : dict
        Dictionary mapping unique atom label to index of `Atom` in `self.asymmetric_unit`.
    anisotropy : dict
        Dictionary mapping unique atom label to displacement anisotropy resolved at that site.
    cif_sgname : str or None
        Space group name obtained by looking up the value of
        `_space_group_name_Hall`,
        `_symmetry_space_group_name_Hall`,
        `_space_group_name_H-M_alt`,
        `_symmetry_space_group_name_H-M`
        items. ``None`` when neither is defined.
    """

    # static data and methods ------------------------------------------------

    # dictionary set of class methods for translating CIF values
    # to Atom attributes

    # static data and methods ------------------------------------------------

    # dictionary set of class methods for translating CIF values
    # to Atom attributes

    _atom_setters = dict.fromkeys(
        (
            "_tr_ignore",
            "_tr_atom_site_label",
            "_tr_atom_site_type_symbol",
            "_tr_atom_site_fract_x",
            "_tr_atom_site_fract_y",
            "_tr_atom_site_fract_z",
            "_tr_atom_site_cartn_x",
            "_tr_atom_site_cartn_y",
            "_tr_atom_site_cartn_z",
            "_tr_atom_site_U_iso_or_equiv",
            "_tr_atom_site_B_iso_or_equiv",
            "_tr_atom_site_adp_type",
            "_tr_atom_site_thermal_displace_type",
            "_tr_atom_site_occupancy",
            "_tr_atom_site_aniso_U_11",
            "_tr_atom_site_aniso_U_22",
            "_tr_atom_site_aniso_U_33",
            "_tr_atom_site_aniso_U_12",
            "_tr_atom_site_aniso_U_13",
            "_tr_atom_site_aniso_U_23",
            "_tr_atom_site_aniso_B_11",
            "_tr_atom_site_aniso_B_22",
            "_tr_atom_site_aniso_B_33",
            "_tr_atom_site_aniso_B_12",
            "_tr_atom_site_aniso_B_13",
            "_tr_atom_site_aniso_B_23",
        )
    )
    # make _atom_setters case insensitive
    for k in list(_atom_setters.keys()):
        _atom_setters[k] = _atom_setters[k.lower()] = k
    del k

    BtoU = 1.0 / (8 * numpy.pi**2)
    """float: Conversion factor from B values to U values."""

    def _tr_ignore(a, value):
        return

    _tr_ignore = staticmethod(_tr_ignore)

    def _tr_atom_site_label(a, value):
        a.label = str(value)
        # set element when not specified by _atom_site_type_symbol
        if not a.element:
            P_cif._tr_atom_site_type_symbol(a, value)

    _tr_atom_site_label = staticmethod(_tr_atom_site_label)

    # 3 regexp groups for nucleon number, atom symbol, and oxidation state
    _psymb = re.compile(r"(\d+-)?([a-zA-Z]+)(\d[+-])?")

    def _tr_atom_site_type_symbol(a, value):
        rx = P_cif._psymb.match(value)
        smbl = rx and rx.group(0) or value
        smbl = str(smbl)
        a.element = smbl[:1].upper() + smbl[1:].lower()

    _tr_atom_site_type_symbol = staticmethod(_tr_atom_site_type_symbol)

    def _tr_atom_site_fract_x(a, value):
        a.xyz[0] = leading_float(value)

    _tr_atom_site_fract_x = staticmethod(_tr_atom_site_fract_x)

    def _tr_atom_site_fract_y(a, value):
        a.xyz[1] = leading_float(value)

    _tr_atom_site_fract_y = staticmethod(_tr_atom_site_fract_y)

    def _tr_atom_site_fract_z(a, value):
        a.xyz[2] = leading_float(value)

    _tr_atom_site_fract_z = staticmethod(_tr_atom_site_fract_z)

    def _tr_atom_site_cartn_x(a, value):
        a.xyz_cartn[0] = leading_float(value)

    _tr_atom_site_cartn_x = staticmethod(_tr_atom_site_cartn_x)

    def _tr_atom_site_cartn_y(a, value):
        a.xyz_cartn[1] = leading_float(value)

    _tr_atom_site_cartn_y = staticmethod(_tr_atom_site_cartn_y)

    def _tr_atom_site_cartn_z(a, value):
        a.xyz_cartn[2] = leading_float(value)

    _tr_atom_site_cartn_z = staticmethod(_tr_atom_site_cartn_z)

    def _tr_atom_site_U_iso_or_equiv(a, value):
        a.Uisoequiv = leading_float(value)

    _tr_atom_site_U_iso_or_equiv = staticmethod(_tr_atom_site_U_iso_or_equiv)

    def _tr_atom_site_B_iso_or_equiv(a, value):
        a.Uisoequiv = P_cif.BtoU * leading_float(value)

    _tr_atom_site_B_iso_or_equiv = staticmethod(_tr_atom_site_B_iso_or_equiv)

    def _tr_atom_site_adp_type(a, value):
        a.anisotropy = value not in ("Uiso", "Biso")

    _tr_atom_site_adp_type = staticmethod(_tr_atom_site_adp_type)
    _tr_atom_site_thermal_displace_type = _tr_atom_site_adp_type

    def _tr_atom_site_occupancy(a, value):
        a.occupancy = leading_float(value, 1.0)

    _tr_atom_site_occupancy = staticmethod(_tr_atom_site_occupancy)

    def _tr_atom_site_aniso_U_11(a, value):
        a.U11 = leading_float(value)

    _tr_atom_site_aniso_U_11 = staticmethod(_tr_atom_site_aniso_U_11)

    def _tr_atom_site_aniso_U_22(a, value):
        a.U22 = leading_float(value)

    _tr_atom_site_aniso_U_22 = staticmethod(_tr_atom_site_aniso_U_22)

    def _tr_atom_site_aniso_U_33(a, value):
        a.U33 = leading_float(value)

    _tr_atom_site_aniso_U_33 = staticmethod(_tr_atom_site_aniso_U_33)

    def _tr_atom_site_aniso_U_12(a, value):
        a.U12 = leading_float(value)

    _tr_atom_site_aniso_U_12 = staticmethod(_tr_atom_site_aniso_U_12)

    def _tr_atom_site_aniso_U_13(a, value):
        a.U13 = leading_float(value)

    _tr_atom_site_aniso_U_13 = staticmethod(_tr_atom_site_aniso_U_13)

    def _tr_atom_site_aniso_U_23(a, value):
        a.U23 = leading_float(value)

    _tr_atom_site_aniso_U_23 = staticmethod(_tr_atom_site_aniso_U_23)

    def _tr_atom_site_aniso_B_11(a, value):
        a.U11 = P_cif.BtoU * leading_float(value)

    _tr_atom_site_aniso_B_11 = staticmethod(_tr_atom_site_aniso_B_11)

    def _tr_atom_site_aniso_B_22(a, value):
        a.U22 = P_cif.BtoU * leading_float(value)

    _tr_atom_site_aniso_B_22 = staticmethod(_tr_atom_site_aniso_B_22)

    def _tr_atom_site_aniso_B_33(a, value):
        a.U33 = P_cif.BtoU * leading_float(value)

    _tr_atom_site_aniso_B_33 = staticmethod(_tr_atom_site_aniso_B_33)

    def _tr_atom_site_aniso_B_12(a, value):
        a.U12 = P_cif.BtoU * leading_float(value)

    _tr_atom_site_aniso_B_12 = staticmethod(_tr_atom_site_aniso_B_12)

    def _tr_atom_site_aniso_B_13(a, value):
        a.U13 = P_cif.BtoU * leading_float(value)

    _tr_atom_site_aniso_B_13 = staticmethod(_tr_atom_site_aniso_B_13)

    def _tr_atom_site_aniso_B_23(a, value):
        a.U23 = P_cif.BtoU * leading_float(value)

    _tr_atom_site_aniso_B_23 = staticmethod(_tr_atom_site_aniso_B_23)

    def _get_atom_setters(cifloop):
        """Static method for finding translators of CifLoop items to data in `Atom` instance.

        Parameters
        ----------
        cifloop : CifLoop
            Instance of `CifLoop`.

        Returns
        -------
        list
            List of setter functions in the order of `cifloop.keys()`.
        """
        rv = []
        for p in cifloop.keys():
            lcname = "_tr" + p.lower()
            fncname = P_cif._atom_setters.get(lcname, "_tr_ignore")
            f = getattr(P_cif, fncname)
            rv.append(f)
        return rv

    _get_atom_setters = staticmethod(_get_atom_setters)

    # normal methods ---------------------------------------------------------

    def __init__(self, eps=None):
        StructureParser.__init__(self)
        self.format = "cif"
        self.ciffile = None
        self.stru = None
        self.spacegroup = None
        self.eps = eps
        self.eau = None
        self.asymmetric_unit = None
        self.labelindex = {}
        self.anisotropy = {}
        self.cif_sgname = None
        pass

    def parse(self, s):
        """Create `Structure` instance from a string in CIF format.

        Parameters
        ----------
        s : str
            A string in CIF format.

        Returns
        -------
        Structure
            `Structure` instance.

        Raises
        ------
        StructureFormatError
            When the data do not constitute a valid CIF format.
        """
        self.ciffile = None
        self.filename = ""
        fp = io.StringIO(s)
        rv = self._parseCifDataSource(fp)
        return rv

    def parseLines(self, lines):
        """Parse list of lines in CIF format.

        Parameters
        ----------
        lines : list
            List of strings stripped of line terminator.

        Returns
        -------
        Structure
            `Structure` instance.

        Raises
        ------
        StructureFormatError
            When the data do not constitute a valid CIF format.
        """
        s = "\n".join(lines) + "\n"
        return self.parse(s)

    def parseFile(self, filename):
        """Create Structure from an existing CIF file.

        Parameters
        ----------
        filename : str
            Path to structure file.

        Returns
        -------
        Structure
            `Structure` instance.

        Raises
        ------
        StructureFormatError
            When the data do not constitute a valid CIF format.
        IOError
            When the file cannot be opened.
        """
        self.ciffile = None
        self.filename = filename
        rv = self._parseCifDataSource(filename)
        # all good here
        return rv

    def _parseCifDataSource(self, datasource):
        """Open and process CIF data from the specified `datasource`.

        Parameters
        ----------
        datasource : str or a file-like object
            This is used as an argument to the `CifFile` class. The `CifFile`
            instance is stored in `ciffile` attribute of this Parser.

        Returns
        -------
        Structure
            The `Structure` object loaded from the specified data source.

        Raises
        ------
        StructureFormatError
            When the data do not constitute a valid CIF format.
        """
        from CifFile import CifFile, StarError

        self.stru = None
        try:
            with _suppressCifParserOutput():
                # Use `grammar` option to digest values with curly-brackets.
                # Ref: https://bitbucket.org/jamesrhester/pycifrw/issues/19
                self.ciffile = CifFile(datasource, grammar="auto")
                for blockname in self.ciffile.keys():
                    self._parseCifBlock(blockname)
                    # stop after reading the first structure
                    if self.stru is not None:
                        break
        except (StarError, ValueError, IndexError) as err:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            emsg = str(err).strip()
            e = StructureFormatError(emsg)
            raise e.with_traceback(exc_traceback)
        return self.stru

    def _parseCifBlock(self, blockname):
        """Translate CIF file block, skip blocks without `_atom_site_label`.
        Updates data members `stru`, `eau`.

        Parameters
        ----------
        blockname : str
            Name of top level block in `self.ciffile`.
        """
        block = self.ciffile[blockname]
        if "_atom_site_label" not in block:
            return
        # here block contains structure, initialize output data
        self.stru = Structure()
        self.labelindex.clear()
        self.anisotropy.clear()
        # execute specialized block parsers
        self._parse_lattice(block)
        self._parse_atom_site_label(block)
        self._parse_atom_site_aniso_label(block)
        self._parse_space_group_symop_operation_xyz(block)
        return

    def _parse_lattice(self, block):
        """Obtain `lattice` parameters from a `CifBlock`.

        This method updates `self.stru.lattic`e.

        Parameters
        ----------
        block : CifBlock
            Instance of CifBlock.
        """
        if "_cell_length_a" not in block:
            return
        # obtain lattice parameters
        try:
            latpars = (
                leading_float(block["_cell_length_a"]),
                leading_float(block["_cell_length_b"]),
                leading_float(block["_cell_length_c"]),
                leading_float(block["_cell_angle_alpha"]),
                leading_float(block["_cell_angle_beta"]),
                leading_float(block["_cell_angle_gamma"]),
            )
        except KeyError as err:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            emsg = str(err)
            e = StructureFormatError(emsg)
            raise e.with_traceback(exc_traceback)
        self.stru.lattice = Lattice(*latpars)
        return

    def _parse_atom_site_label(self, block):
        """Obtain atoms in asymmetric unit from a `CifBlock`.

        This method inserts `Atom` instances to `self.stru` and
        updates `labelindex` dictionary.

        Parameters
        ----------
        block : CifBlock
            Instance of `CifBlock`.
        """
        # process _atom_site_label
        atom_site_loop = block.GetLoop("_atom_site_label")
        does_adp_type = (
            "_atom_site_adp_type" in atom_site_loop or "_atom_site_thermal_displace_type" in atom_site_loop
        )
        # get a list of setters for atom_site values
        prop_setters = P_cif._get_atom_setters(atom_site_loop)
        # index of the _atom_site_label item for the labelindex dictionary
        ilb = atom_site_loop.keys().index("_atom_site_label")
        # loop through the values and pass them to the setters
        sitedatalist = zip(*atom_site_loop.values())
        for values in sitedatalist:
            curlabel = values[ilb]
            # skip entries that have invalid label
            if curlabel == "?":
                continue
            self.labelindex[curlabel] = len(self.stru)
            self.stru.addNewAtom()
            a = self.stru.getLastAtom()
            for fset, val in zip(prop_setters, values):
                fset(a, val)
            if does_adp_type:
                self.anisotropy[curlabel] = a.anisotropy
        return

    def _parse_atom_site_aniso_label(self, block):
        """Obtain value of anisotropic thermal displacements from a `CifBlock`.

        This method updates `U` members of `Atom` instances in `self.stru`.
        The `labelindex` dictionary has to be defined beforehand.

        Parameters
        ----------
        block : CifBlock
            Instance of `CifBlock`.
        """
        if "_atom_site_aniso_label" not in block:
            return
        # something to do here:
        adp_loop = block.GetLoop("_atom_site_aniso_label")
        # index of the _atom_site_label column
        ilb = adp_loop.keys().index("_atom_site_aniso_label")
        # get a list of setters for this loop
        prop_setters = P_cif._get_atom_setters(adp_loop)
        sitedatalist = zip(*adp_loop.values())
        for values in sitedatalist:
            lb = values[ilb]
            if lb == "?":
                break
            idx = self.labelindex[lb]
            a = self.stru[idx]
            if lb not in self.anisotropy:
                a.anisotropy = True
                self.anisotropy[lb] = True
            for fset, val in zip(prop_setters, values):
                fset(a, val)
        return

    def _parse_space_group_symop_operation_xyz(self, block):
        """Process symmetry operations from a CifBlock.

        The method updates `spacegroup` and `eau` data according to symmetry
        operations defined in `_space_group_symop_operation_xyz` or
        `_symmetry_equiv_pos_as_xyz` items in `CifBlock`.

        Parameters
        ----------
        block : CifBlock
            Instance of `CifBlock`.
        """
        from diffpy.structure.spacegroups import FindSpaceGroup, GetSpaceGroup, IsSpaceGroupIdentifier, SpaceGroup

        self.asymmetric_unit = list(self.stru)
        sym_synonyms = ("_space_group_symop_operation_xyz", "_symmetry_equiv_pos_as_xyz")
        sym_loop_name = [n for n in sym_synonyms if n in block]
        # recover explicit list of symmetry operations
        symop_list = []
        if sym_loop_name:
            # sym_loop exists here and we know its cif name
            sym_loop_name = sym_loop_name[0]
            sym_loop = block.GetLoop(sym_loop_name)
            for eqxyz in sym_loop[sym_loop_name]:
                opcif = getSymOp(eqxyz)
                symop_list.append(opcif)
        # determine space group number
        sg_nameHall = block.get("_space_group_name_Hall", "") or block.get("_symmetry_space_group_name_Hall", "")
        sg_nameHM = (
            block.get("_space_group_name_H-M_alt", "")
            or block.get("_space_group_name_H-M_ref", "")
            or block.get("_symmetry_space_group_name_H-M", "")
        )
        self.cif_sgname = sg_nameHall or sg_nameHM or None
        sgid = block.get("_space_group_IT_number", "") or block.get("_symmetry_Int_Tables_number", "") or sg_nameHM
        self.spacegroup = None
        # try to reuse existing space group from symmetry operations
        if symop_list:
            try:
                self.spacegroup = FindSpaceGroup(symop_list)
            except ValueError:
                pass
        # otherwise lookup the space group from its identifier
        if self.spacegroup is None and sgid and IsSpaceGroupIdentifier(sgid):
            self.spacegroup = GetSpaceGroup(sgid)
        # define new spacegroup when symmetry operations were listed, but
        # there is no match to an existing definition
        if symop_list and self.spacegroup is None:
            new_short_name = "CIF " + (sg_nameHall or "data")
            new_crystal_system = (
                block.get("_space_group_crystal_system") or block.get("_symmetry_cell_setting") or "TRICLINIC"
            ).upper()
            self.spacegroup = SpaceGroup(
                short_name=new_short_name, crystal_system=new_crystal_system, symop_list=symop_list
            )
        if self.spacegroup is None:
            emsg = "CIF file has unknown space group identifier {!r}."
            raise StructureFormatError(emsg.format(sgid))
        self._expandAsymmetricUnit(block)
        return

    def _expandAsymmetricUnit(self, block):
        """Perform symmetry expansion of `self.stru` using `self.spacegroup`.

        This method updates data in `stru` and `eau`.

        Parameters
        ----------
        block : CifBlock
            The top-level block containing crystal structure data.
        """
        from diffpy.structure.symmetryutilities import ExpandAsymmetricUnit

        corepos = [a.xyz for a in self.stru]
        coreUijs = [a.U for a in self.stru]
        self.eau = ExpandAsymmetricUnit(self.spacegroup, corepos, coreUijs, eps=self.eps)
        # setup anisotropy according to symmetry requirements
        # unless it was already explicitly set
        for ca, uisotropy in zip(self.stru, self.eau.Uisotropy):
            if ca.label not in self.anisotropy:
                ca.anisotropy = not uisotropy
                self.anisotropy[ca.label] = ca.anisotropy
        # build a nested list of new atoms:
        newatoms = []
        for i, ca in enumerate(self.stru):
            eca = []  # expanded core atom
            for j in range(self.eau.multiplicity[i]):
                a = Atom(ca)
                a.xyz = self.eau.expandedpos[i][j]
                if j > 0:
                    a.label += "_" + str(j + 1)
                if a.anisotropy:
                    a.U = self.eau.expandedUijs[i][j]
                eca.append(a)
            newatoms.append(eca)
        # insert new atoms where they belong
        self.stru[:] = sum(newatoms, [])
        return

    # conversion to CIF ------------------------------------------------------

    def toLines(self, stru):
        """Convert `Structure` to a list of lines in basic CIF format.

        Parameters
        ----------
        stru : Structure
            The structure to be converted.

        Returns
        -------
        list
            List of lines in basic CIF format.
        """
        import time

        lines = []
        # may be replaced with filtered Structure.title
        # for now, we can add the title as a comment
        if stru.title.strip() != "":
            title_lines = stru.title.split("\n")
            lines.extend(["# " + line.strip() for line in title_lines])
            lines.append("")
        lines.append("data_3D")
        iso_date = "%04i-%02i-%02i" % time.gmtime()[:3]
        lines.extend(
            [
                "%-31s %s" % ("_audit_creation_date", iso_date),
                "%-31s %s" % ("_audit_creation_method", "P_cif.py"),
                "",
                "%-31s %s" % ("_symmetry_space_group_name_H-M", "'P1'"),
                "%-31s %s" % ("_symmetry_Int_Tables_number", "1"),
                "%-31s %s" % ("_symmetry_cell_setting", "triclinic"),
                "",
            ]
        )
        # there should be no need to specify equivalent positions for P1
        # _symmetry_equiv_posi_as_xyz x,y,z
        lines.extend(
            [
                "%-31s %.6g" % ("_cell_length_a", stru.lattice.a),
                "%-31s %.6g" % ("_cell_length_b", stru.lattice.b),
                "%-31s %.6g" % ("_cell_length_c", stru.lattice.c),
                "%-31s %.6g" % ("_cell_angle_alpha", stru.lattice.alpha),
                "%-31s %.6g" % ("_cell_angle_beta", stru.lattice.beta),
                "%-31s %.6g" % ("_cell_angle_gamma", stru.lattice.gamma),
                "",
            ]
        )
        # build a list of site labels and adp (displacement factor) types
        element_count = {}
        a_site_label = []
        a_adp_type = []
        for a in stru:
            cnt = element_count[a.element] = element_count.get(a.element, 0) + 1
            a_site_label.append("%s%i" % (a.element, cnt))
            if numpy.all(a.U == a.U[0, 0] * numpy.identity(3)):
                a_adp_type.append("Uiso")
            else:
                a_adp_type.append("Uani")
        # list all atoms
        lines.extend(
            [
                "loop_",
                "  _atom_site_label",
                "  _atom_site_type_symbol",
                "  _atom_site_fract_x",
                "  _atom_site_fract_y",
                "  _atom_site_fract_z",
                "  _atom_site_U_iso_or_equiv",
                "  _atom_site_adp_type",
                "  _atom_site_occupancy",
            ]
        )
        for i in range(len(stru)):
            a = stru[i]
            line = "  %-5s %-3s %11.6f %11.6f %11.6f %11.6f %-5s %.4f" % (
                a_site_label[i],
                a.element,
                a.xyz[0],
                a.xyz[1],
                a.xyz[2],
                a.Uisoequiv,
                a_adp_type[i],
                a.occupancy,
            )
            lines.append(line)
        # find anisotropic atoms
        idx_aniso = [i for i in range(len(stru)) if a_adp_type[i] != "Uiso"]
        if idx_aniso != []:
            lines.extend(
                [
                    "loop_",
                    "  _atom_site_aniso_label",
                    "  _atom_site_aniso_U_11",
                    "  _atom_site_aniso_U_22",
                    "  _atom_site_aniso_U_33",
                    "  _atom_site_aniso_U_12",
                    "  _atom_site_aniso_U_13",
                    "  _atom_site_aniso_U_23",
                ]
            )
            for i in idx_aniso:
                a = stru[i]
                line = "  %-5s %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f" % (
                    a_site_label[i],
                    a.U[0, 0],
                    a.U[1, 1],
                    a.U[2, 2],
                    a.U[0, 1],
                    a.U[0, 2],
                    a.U[1, 2],
                )
                lines.append(line)
        return lines


# End of class P_cif

# Routines -------------------------------------------------------------------

# constant regular expression for leading_float()
rx_float = re.compile(r"[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?")


def leading_float(s, d=0.0):
    """Extract the first float from a string and ignore trailing characters.

    Useful for extracting values from "value(std)" syntax.

    Parameters
    ----------
    s : str
        The string to be scanned for floating point value.
    d : float, Optional
        The default value when `s` is "." or "?", which in CIF
        format stands for inapplicable and unknown, respectively.

    Returns
    -------
    float
        The extracted floating point value.

    Raises
    ------
    ValueError
        When string does not start with a float.
    """
    sbare = s.strip()
    mx = rx_float.match(sbare)
    if mx:
        rv = float(mx.group())
    elif sbare == "." or sbare == "?":
        # CIF files may contain "." or "?" for unknown values
        rv = d
    else:
        rv = float(sbare)
    return rv


# helper dictionary for getSymOp()
symvec = {
    "x": numpy.array([1, 0, 0], dtype=float),
    "y": numpy.array([0, 1, 0], dtype=float),
    "z": numpy.array([0, 0, 1], dtype=float),
    "-x": numpy.array([-1, 0, 0], dtype=float),
    "-y": numpy.array([0, -1, 0], dtype=float),
    "-z": numpy.array([0, 0, -1], dtype=float),
}
symvec["+x"] = symvec["x"]
symvec["+y"] = symvec["y"]
symvec["+z"] = symvec["z"]


def getSymOp(s):
    """Create `SpaceGroups.SymOp` instance from a string.

    Parameters
    ----------
    s : str
        Formula for equivalent coordinates, for example ``'x,1/2-y,1/2+z'``.

    Returns
    -------
    SymOp
        Instance of `SymOp`.
    """
    from diffpy.structure.spacegroups import SymOp

    snoblanks = s.replace(" ", "")
    eqlist = snoblanks.split(",")
    R = numpy.zeros((3, 3), dtype=float)
    t = numpy.zeros(3, dtype=float)
    for i in (0, 1, 2):
        eqparts = re.split("(?i)([+-]?[xyz])", eqlist[i])
        for Rpart in eqparts[1::2]:
            R[i, :] += symvec[Rpart.lower()]
        for tpart in eqparts[::2]:
            t[i] += eval("1.0*%s+0" % tpart)
    t -= numpy.floor(t)
    rv = SymOp(R, t)
    return rv


def getParser(eps=None):
    """Return new `parser` object for CIF format.

    Parameters
    ----------
    eps : float, Optional
        fractional coordinates cutoff for duplicate positions.
        When ``None`` use the default for `ExpandAsymmetricUnit`: ``1.0e-5``.

    Returns
    -------
    P_cif
        Instance of `P_cif`.
    """
    return P_cif(eps=eps)


# Local Helpers --------------------------------------------------------------


@contextmanager
def _suppressCifParserOutput():
    """Context manager which suppresses diagnostic messages from CIF parser."""
    from CifFile import yapps3_compiled_rt

    print_error = yapps3_compiled_rt.print_error
    # replace the print_error function with no-operation
    yapps3_compiled_rt.print_error = lambda *a, **kw: None
    try:
        yield print_error
    finally:
        yapps3_compiled_rt.print_error = print_error
    pass

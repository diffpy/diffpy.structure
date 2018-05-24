#!/usr/bin/env python
##############################################################################
#
# diffpy.Structure  by DANSE Diffraction group
#                   Simon J. L. Billinge
#                   (c) 2011 Trustees of the Columbia University
#                   in the City of New York.  All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE_DANSE.txt for license information.
#
##############################################################################


'''Space group classes and definitions from mmLib and sgtbx.
'''

from diffpy.Structure.spacegroupmod import *
from diffpy.Structure.mmlibspacegroups import *
from diffpy.Structure.sgtbxspacegroups import *


# all spacegroup definitions
SpaceGroupList = mmLibSpaceGroupList + sgtbxSpaceGroupList


def GetSpaceGroup(sgid):
    """Returns the SpaceGroup instance for the given identifier.

    sgid -- space group symbol, either short_name or pdb_name,
            whatever it means in mmlib.  Can be also an integer.

    Return space group instance.
    Raise ValueError when not found.
    """
    if not _sg_lookup_table:
        _buildSGLookupTable()
    if sgid in _sg_lookup_table:
        return _sg_lookup_table[sgid]
    # Try different versions of sgid, first make sure it is a string
    emsg = "Unknown space group identifier %r" % sgid
    if not isinstance(sgid, basestring):
        raise ValueError(emsg)
    # short name case adjusted
    sgkey = sgid.strip()
    sgkey = sgkey[:1].upper() + sgkey[1:].lower()
    if sgkey in _sg_lookup_table:
        return _sg_lookup_table[sgkey]
    # long name all upper case
    sgkey = sgid.strip().upper()
    if sgkey in _sg_lookup_table:
        return _sg_lookup_table[sgkey]
    # try to remove any blanks
    sgkey = sgid.replace(' ', '')
    if sgkey in _sg_lookup_table:
        return _sg_lookup_table[sgkey]
    # nothing worked, sgid is unknown identifier
    raise ValueError(emsg)


def IsSpaceGroupIdentifier(sgid):
    """Check if identifier can be used as an argument to GetSpaceGroup.

    Return bool.
    """
    try:
        GetSpaceGroup(sgid)
        rv = True
    except ValueError:
        rv = False
    return rv


def _buildSGLookupTable():
    """Rebuild space group lookup table from the SpaceGroupList data.

    This routine updates the global _sg_lookup_table dictionary.
    No return value.
    """
    _sg_lookup_table.clear()
    for sg in SpaceGroupList:
        _sg_lookup_table.setdefault(sg.number, sg)
        _sg_lookup_table.setdefault(str(sg.number), sg)
        _sg_lookup_table.setdefault(sg.short_name, sg)
        _sg_lookup_table.setdefault(sg.pdb_name, sg)
        _sg_lookup_table.setdefault(sg.alt_name, sg)
    # extra aliases obtained from matching code in
    # cctbx::sgtbx::symbols::find_main_symbol_dict_entry
    alias_hmname = [
        ('Pm3', 'P m -3'),
        ('Pn3', 'P n -3'),
        ('Fm3', 'F m -3'),
        ('Fd3', 'F d -3'),
        ('Im3', 'I m -3'),
        ('Pa3', 'P a -3'),
        ('Ia3', 'I a -3'),
        ('Pm3m', 'P m -3 m'),
        ('Pn3n', 'P n -3 n'),
        ('Pm3n', 'P m -3 n'),
        ('Pn3m', 'P n -3 m'),
        ('Fm3m', 'F m -3 m'),
        ('Fm3c', 'F m -3 c'),
        ('Fd3m', 'F d -3 m'),
        ('Fd3c', 'F d -3 c'),
        ('Im3m', 'I m -3 m'),
        ('Ia3d', 'I a -3 d'),
    ]
    for a, hm in alias_hmname:
        hmbare = hm.replace(' ', '')
        _sg_lookup_table.setdefault(a, _sg_lookup_table[hmbare])
    # make sure None does not sneak into the dictionary
    if None in _sg_lookup_table:
        del _sg_lookup_table[None]
    return
_sg_lookup_table = {}


# End of file

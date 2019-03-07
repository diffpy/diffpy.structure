#!/usr/bin/env python

'''Quick and extremely dirty script for generating code for SpaceGroup that
are defined in cctbx, but not in mmLib.  It was used to generate module
sgtbxspacegroups.

This is a utility script that should not be included with code distribution.

Not to be included with code distributions.
'''


import re
import math
import numpy
from diffpy.structure.spacegroups import SpaceGroup, SymOp
from diffpy.structure.spacegroups import mmLibSpaceGroupList
from diffpy.structure.spacegroups import IsSpaceGroupIdentifier
from cctbx import sgtbx

def tupleToSGArray(tpl):
    if not _rtarrays:
        import diffpy.structure.SpaceGroups as sgmod
        for n in dir(sgmod):
            if not n.startswith('Rot_') and not n.startswith('Tr_'): continue
            a = getattr(sgmod, n)
            t = tuple(a.flatten())
            _rtarrays[t] = a
    if len(tpl) == 3:
        tpl = tuple([(x - math.floor(x)) for x in tpl])
        if tpl not in _rtarrays:
            _rtarrays[tpl] = numpy.array(tpl, dtype=float)
    return _rtarrays[tpl]
_rtarrays = {}


def mmSpaceGroupFromSymbol(symbol):
    """Construct SpaceGroup instance from a string symbol using sgtbx data.
    """
    sginfo = sgtbx.space_group_info(symbol)
    symop_list = []
    symop_list = getSymOpList(sginfo.group())
    sgtype = sginfo.type()
    uhm = sgtype.lookup_symbol()
    sgsmbls = sgtbx.space_group_symbols(uhm)
    kw = {}
    kw['number'] = sgtype.number()
    kw['num_sym_equiv'] = len(symop_list)
    kw['num_primitive_sym_equiv'] = countUniqueRotations(symop_list)
    kw['short_name'] = sgsmbls.hermann_mauguin().replace(' ', '')
    pgt = sgsmbls.point_group_type()
    pgn = "PG" + re.sub('-(\d)', '\\1bar', pgt)
    kw['point_group_name'] = pgn
    kw['crystal_system'] = sgsmbls.crystal_system().upper()
    kw['pdb_name'] = sgsmbls.hermann_mauguin()
    kw['symop_list'] = symop_list
    mmsg = SpaceGroup(**kw)
    return mmsg


def adjustMMSpaceGroupNumber(mmsg):
    sg0 = [x for x in mmLibSpaceGroupList if x.number == mmsg.number]
    if sg0 and cmpSpaceGroups(sg0[0], mmsg):
        return
    while mmsg.number in sgnumbers:
        mmsg.number += 1000
    sgnumbers.append(mmsg.number)


def getSymOpList(grp):
    symop_list = []
    for op in grp:
        r_sgtbx = op.r().as_double()
        t_sgtbx = op.t().as_double()
        R = tupleToSGArray(r_sgtbx)
        t = tupleToSGArray(t_sgtbx)
        symop_list.append(SymOp(R, t))
    return symop_list


def countUniqueRotations(symop_list):
    unique_rotations = set()
    for op in symop_list:
        tpl = tuple(op.R.flatten())
        unique_rotations.add(tpl)
    return len(unique_rotations)


def cmpSpaceGroups(sg0, sg1):
    if sg0 is sg1:  return True
    s0 = hashMMSpaceGroup(sg0)
    s1 = hashMMSpaceGroup(sg1)
    return s0 == s1


def findEquivalentMMSpaceGroup(grp):
    if not _equivmmsg:
        for sgn in mmLibSpaceGroupList:
            ssgn = hashMMSpaceGroup(sgn)
            _equivmmsg.setdefault(ssgn, sgn)
    ssg = hashSgtbxGroup(grp)
    return _equivmmsg.get(ssg)
_equivmmsg = {}


def findEquivalentSgtbxSpaceGroup(sgmm):
    if not _equivsgtbx:
        for smbls in sgtbx.space_group_symbol_iterator():
            uhm = smbls.universal_hermann_mauguin()
            grp = sgtbx.space_group_info(uhm).group()
            hgrp = hashSgtbxGroup(grp)
            _equivsgtbx.setdefault(hgrp, grp)
    hgmm = hashMMSpaceGroup(sgmm)
    return _equivsgtbx.get(hgmm)
_equivsgtbx = {}


def hashMMSpaceGroup(sg):
    lines = [str(sg.number % 1000)] + sorted(map(str, sg.iter_symops()))
    s = '\n'.join(lines)
    return s


def hashSgtbxGroup(grp):
    n = grp.type().number()
    lines = [str(n)] + sorted(map(str, getSymOpList(grp)))
    s = '\n'.join(lines)
    return s

sgnumbers = [sg.number for sg in mmLibSpaceGroupList]

_SGsrc = '''\
sg%(number)i = SpaceGroup(
    number = %(number)i,
    num_sym_equiv = %(num_sym_equiv)i,
    num_primitive_sym_equiv = %(num_sym_equiv)i,
    short_name = %(short_name)r,
    point_group_name = %(point_group_name)r,
    crystal_system = %(crystal_system)r,
    pdb_name = %(pdb_name)r,
    symop_list = [
        @SYMOPS@
    ]
)
'''

def SGCode(mmsg):
    src0 = _SGsrc % mmsg.__dict__
    src1 = src0.replace('@SYMOPS@', SymOpsCode(mmsg))
    return src1


def SymOpsCode(mmsg):
    lst = ["%8s%s," % ('', SymOpCode(op)) for op in mmsg.iter_symops()]
    src = '\n'.join(lst).strip()
    return src


def SymOpCode(op):
    if not _rtnames:
        import diffpy.structure.SpaceGroups as sgmod
        for n in dir(sgmod):
            if not n.startswith('Rot_') and not n.startswith('Tr_'): continue
            a = getattr(sgmod, n)
            at = tuple(a.flatten())
            _rtnames[at] = 'sgmod.' + n
    nR = _rtnames[tuple(op.R.flatten())]
    nt = _rtnames[tuple(op.t)]
    src = 'SymOp(%s, %s)' % (nR, nt)
    return src
_rtnames = {}


def main():
    duplicates = set()
    for smbls in sgtbx.space_group_symbol_iterator():
        uhm = smbls.universal_hermann_mauguin()
        grp = sgtbx.space_group_info(uhm).group()
        if findEquivalentMMSpaceGroup(grp): continue
        shn = smbls.hermann_mauguin().replace(' ', '')
        if IsSpaceGroupIdentifier(shn): continue
        sg = mmSpaceGroupFromSymbol(uhm)
        hsg = hashMMSpaceGroup(sg)
        if hsg in duplicates: continue
        adjustMMSpaceGroupNumber(sg)
        duplicates.add(hsg)
        print(SGCode(sg))
    return

if __name__ == '__main__':
    main()

#!/usr/bin/env python
##############################################################################
#
# diffpy.structure  Complex Modeling Initiative
#                   (c) 2016 Brookhaven Science Associates,
#                   Brookhaven National Laboratory.
#                   All rights reserved.
#
# File coded by:    Pavol Juhas
#
# See AUTHORS.txt for a list of people who contributed.
# See LICENSE.txt for license information.
#
##############################################################################
"""Test for Structure utilities."""
import pytest

from diffpy.structure.utils import atom_bare_symbol, atomBareSymbol


def test_atomBareSymbol():
    assert atomBareSymbol("Cl-") == "Cl"
    assert atomBareSymbol("Ca2+") == "Ca"
    assert atomBareSymbol("12-C") == "C"


@pytest.mark.parametrize(
    "symbol, expected",
    [
        ("Cl-", "Cl"),
        ("Ca2+", "Ca"),
        ("12-C", "C"),
    ],
)
def test_atom_bare_symbol(symbol, expected):
    actual = atom_bare_symbol(symbol)
    assert actual == expected

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
# See AUTHORS.rst for a list of people who contributed.
# See LICENSE.rst for license information.
#
##############################################################################
"""Tests for the expansion module utilities."""
import pytest

from diffpy.structure.atom import Atom
from diffpy.structure.expansion.makeellipsoid import make_ellipsoid, make_sphere, makeEllipsoid, makeSphere
from diffpy.structure.expansion.shapeutils import find_center, findCenter
from diffpy.structure.lattice import Lattice
from diffpy.structure.structure import Structure


def test_findCenter():
    # C1: We have single atom, expect to return the index of the single atom.
    structure_1 = Structure([Atom("Ni", [0.8, 1.2, 0.9])], lattice=Lattice())
    expected = 0
    actual = findCenter(structure_1)
    assert actual == expected

    # C2: We have multiple atoms.
    # Expect to find the index of the atom which has the closest distance to [0.5, 0.5, 0.5].
    # In this case it corresponds to atom "C".
    structure_2 = Structure(
        [Atom("Ni", [0.8, 1.2, 0.9]), Atom("C", [0.1, 0.2, 0.3])],
        lattice=Lattice(),
    )
    actual = findCenter(structure_2)
    expected = 1
    assert actual == expected


def test_makeEllipsoid_and_makeSphere():
    structure = Structure(
        [
            Atom("Ni", [0.0, 0.0, 0.0]),
            Atom("C", [0.5, 0.5, 0.5]),
            Atom("O", [1.0, 0.0, 0.0]),
            Atom("H", [1.1, 0.0, 0.0]),
        ],
        lattice=Lattice(1, 1, 1, 90, 90, 90),
    )
    # C1: set primary, secondary, and polar radius the same in makeEllipsoid, expect to be the same as makeSphere
    ellipsoid_1 = makeEllipsoid(structure, 1, 1, 1)
    sphere = makeSphere(structure, 1)
    assert [atom.element for atom in ellipsoid_1] == [atom.element for atom in sphere]
    assert [tuple(atom.xyz) for atom in ellipsoid_1] == [tuple(atom.xyz) for atom in sphere]
    # C2: set the radius to be 0.5, expect to exclude the atom that is too far from the center of ellipsoid.
    ellipsoid_2 = makeEllipsoid(structure, 0.5)
    actual = [(atom.element, tuple(atom.xyz)) for atom in ellipsoid_2]
    expected = [("C", (0.5, 0.5, 0.5))]
    assert actual == expected


@pytest.mark.parametrize(
    "atom_params, expected",
    [
        pytest.param([Atom("Ni", [0.8, 1.2, 0.9])], 0),
        pytest.param([Atom("Ni", [0.8, 1.2, 0.9]), Atom("C", [0.1, 0.2, 0.3])], 1),
    ],
)
def test_find_center(atom_params, expected):
    structure = Structure(atom_params, lattice=Lattice())
    actual = find_center(structure)
    assert actual == expected


@pytest.mark.parametrize(
    ("ellipsoid_args", "sphere_radius"),
    [
        ((1, 1, 1), 1),
        ((0.8, 0.8, 0.8), 0.8),
        ((0.5, 0.5, 0.5), 0.5),
    ],
)
def test_make_ellipsoid_equiv_to_make_sphere(ellipsoid_args, sphere_radius):
    structure = Structure(
        [
            Atom("Ni", [0.0, 0.0, 0.0]),
            Atom("C", [0.5, 0.5, 0.5]),
            Atom("O", [1.0, 0.0, 0.0]),
            Atom("H", [1.1, 0.0, 0.0]),
        ],
        lattice=Lattice(1, 1, 1, 90, 90, 90),
    )

    actual = [(atom.element, tuple(atom.xyz)) for atom in make_ellipsoid(structure, *ellipsoid_args)]
    expected = [(atom.element, tuple(atom.xyz)) for atom in make_sphere(structure, sphere_radius)]

    assert actual == expected


@pytest.mark.parametrize(
    ("ellipsoid_args", "expected"),
    [
        ((0.5,), [("C", (0.5, 0.5, 0.5))]),
        ((0.4,), [("C", (0.5, 0.5, 0.5))]),
    ],
)
def test_make_ellipsoid(ellipsoid_args, expected):
    structure = Structure(
        [
            Atom("Ni", [0.0, 0.0, 0.0]),
            Atom("C", [0.5, 0.5, 0.5]),
            Atom("O", [1.0, 0.0, 0.0]),
            Atom("H", [1.1, 0.0, 0.0]),
        ],
        lattice=Lattice(1, 1, 1, 90, 90, 90),
    )
    actual = [(atom.element, tuple(atom.xyz)) for atom in make_ellipsoid(structure, *ellipsoid_args)]
    assert actual == expected

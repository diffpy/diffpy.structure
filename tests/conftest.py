import json
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms

from diffpy.structure import Atom, Lattice, Structure


@pytest.fixture
def user_filesystem(tmp_path):
    base_dir = Path(tmp_path)
    home_dir = base_dir / "home_dir"
    home_dir.mkdir(parents=True, exist_ok=True)
    cwd_dir = base_dir / "cwd_dir"
    cwd_dir.mkdir(parents=True, exist_ok=True)

    home_config_data = {"username": "home_username", "email": "home@email.com"}
    with open(home_dir / "diffpyconfig.json", "w") as f:
        json.dump(home_config_data, f)

    yield tmp_path


@pytest.fixture
def datafile():
    """Fixture to dynamically load any test file."""

    def _load(filename):
        return "tests/testdata/" + filename

    return _load


@pytest.fixture
def build_ase_atom_object():
    """Helper function to build an ASE.Atoms object for testing."""
    a = 5.409
    frac_coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.5],
            [0.25, 0.25, 0.25],
            [0.75, 0.75, 0.75],
        ]
    )
    cart_coords = frac_coords * a
    symbols = ["Zn", "Zn", "S", "S"]
    ase_zb = Atoms(symbols=symbols, positions=cart_coords, cell=[[a, 0, 0], [0, a, 0], [0, 0, a]], pbc=True)
    return ase_zb


@pytest.fixture
def build_diffpy_structure_object():
    """Helper function to build a diffpy.structure.Structure object for
    testing."""
    a = 5.409
    frac_coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.5],
            [0.25, 0.25, 0.25],
            [0.75, 0.75, 0.75],
        ]
    )
    lattice = Lattice(base=[[a, 0, 0], [0, a, 0], [0, 0, a]])
    atoms = [
        Atom("Zn", frac_coords[0]),
        Atom("Zn", frac_coords[1]),
        Atom("S", frac_coords[2]),
        Atom("S", frac_coords[3]),
    ]
    diffpy_zb = Structure(atoms=atoms, lattice=lattice)
    return diffpy_zb

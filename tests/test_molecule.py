from PyVOTCA.molecule import Molecule
from .utils import PATH_TEST
from pathlib import Path
import numpy as np


def test_molecule_IO(tmp_path: Path):
    """Check molecule methods."""
    mol = Molecule()

    # check reading method
    mol.readXYZfile(PATH_TEST / "ethylene.xyz")
    assert len(mol.elements) == 6

    # Check writing method
    file_name = tmp_path / "test.xyz"
    mol.writeXYZfile(file_name)
    assert file_name.exists()

def test_atom_add():
    """"Check that atoms are added properly."""
    mol = Molecule()
    mol.add_atom("O", 1.2,1.44,5.32)
    assert mol.elements[0] == "O"
    assert np.allclose(mol.coordinates[0], [1.2, 1.44, 5.32])

def test_total_energies():
    """Check that total energies are reported correctly."""
    mol = Molecule()
    mol.readORB('example.orb')
    dft_en = mol.getTotalEnergy(kind='dft_tot')
    dft_en_ref = -113.21956697102704
    assert np.isclose(dft_en,dft_en_ref)
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
    mol.add_atom("O", 1.2, 1.44, 5.32)
    assert mol.elements[0] == "O"
    assert np.allclose(mol.coordinates[0], [1.2, 1.44, 5.32])


def test_total_energies():
    """Check that total energies are reported correctly."""
    mol = Molecule()
    mol.readORB(PATH_TEST / "example.orb")
    # total DFT enery
    dft_en = mol.getDFTEnergy()
    dft_en_ref = -113.21956697102704
    assert np.isclose(dft_en, dft_en_ref)
    # KS total energies
    ks_homo_en_ref = -112.81958947
    ks_homo_en = mol.getKSTotalEnergy(6)
    assert np.isclose(ks_homo_en, ks_homo_en_ref)
    ks_lumo_en_ref = -113.26985321
    ks_lumo_en = mol.getKSTotalEnergy(7)
    assert np.isclose(ks_lumo_en, ks_lumo_en_ref)
    # QPpert total energies
    qpp_homo_en_ref = -112.69993863
    qpp_homo_en = mol.getQPTotalEnergy(6)
    assert np.isclose(qpp_homo_en, qpp_homo_en_ref)
    qpp_lumo_en_ref = -113.13870327
    qpp_lumo_en = mol.getQPTotalEnergy(level=7)
    assert np.isclose(qpp_lumo_en, qpp_lumo_en_ref)
    # QPdiag total energies
    qpd_homo_en_ref = -112.7000959
    qpd_homo_en = mol.getQPdiagTotalEnergy(6)
    assert np.isclose(qpd_homo_en, qpd_homo_en_ref)
    qpd_lumo_en_ref = -113.14263116
    qpd_lumo_en = mol.getQPdiagTotalEnergy(level=7)
    assert np.isclose(qpd_lumo_en, qpd_lumo_en_ref)


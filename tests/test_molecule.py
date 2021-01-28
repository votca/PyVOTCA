from pyvotca.molecule import Molecule
from .utils import PATH_TEST
from pathlib import Path
import numpy as np
import unittest


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
    # BSE Singlet total energies
    singlet_en_ref = -112.93707785
    singlet_en = mol.getBSEsingletTotalEnergy(0)
    assert np.isclose(singlet_en, singlet_en_ref)
    # BSE triplet total energies
    triplet_en_ref = -113.02898113
    triplet_en = mol.getBSEtripletTotalEnergy(0)
    assert np.isclose(triplet_en, triplet_en_ref)
    # BSE dynamic Singlet total energies
    singlet_dyn_en_ref = -112.94138747
    singlet_dyn_en = mol.getBSEsingletDynamicTotalEnergy(0)
    assert np.isclose(singlet_dyn_en, singlet_dyn_en_ref)
    # BSE dynamic triplet total energies
    triplet_dyn_en_ref = -113.031518836
    triplet_dyn_en = mol.getBSEtripletDynamicTotalEnergy(0)
    assert np.isclose(triplet_dyn_en, triplet_dyn_en_ref)
    # requesting unavailable KS level
    assert(np.isclose(mol.getKSTotalEnergy(1000), 0.0))
    # requesting unavailable QP level
    assert(np.isclose(mol.getQPTotalEnergy(1000), 0.0))
    # requesting unavailable QP diag level
    assert(np.isclose(mol.getQPdiagTotalEnergy(1000), 0.0))
    # requesting unavailable BSE singlet level
    assert(np.isclose(mol.getBSEsingletTotalEnergy(1000), 0.0))
    # requesting unavailable BSE triplet level
    assert(np.isclose(mol.getBSEtripletTotalEnergy(1000), 0.0))
    # requesting unavailable BSE singlet dynamic level
    assert(np.isclose(mol.getBSEsingletDynamicTotalEnergy(1000), 0.0))
    # requesting unavailable BSE triplet dynamic level
    assert(np.isclose(mol.getBSEtripletDynamicTotalEnergy(1000), 0.0))


def test_QP_corrections():
    """Check that QP corrections are reported correctly."""
    mol = Molecule()
    mol.readORB(PATH_TEST / "example.orb")
    qp_corr_ref = np.array([-0.66958499, -0.58588884, -0.10798883, -0.12099109, -0.11734361, -0.11734361, -0.11965084, 0.13114994, 0.13114994,
                            0.0762104, 0.07464675, 0.0875364, 0.0875364, 0.09781442, 0.10893873, 0.10893873, 0.14590749, 0.12154184, 0.1215417, 0.10959073])
    qp_corr = mol.getQPcorrections()
    assert(np.allclose(qp_corr, qp_corr_ref))


def test_oscillatorStrengths():
    """Check that oscillator strengths  are reported correctly."""
    mol = Molecule()
    mol.readORB(PATH_TEST / "example.orb")
    e_ref = np.array(
        [0.28248912, 0.28248905, 0.2903493,  0.31353545, 0.31353545])
    osc_ref = np.array([3.93892214e-02, 3.93892190e-02, 7.26550288e-10, 4.52526432e-16,
                        4.52526432e-16])
    e, osc = mol.getOscillatorStrengths()
    assert(np.allclose(e, e_ref))
    assert(np.allclose(osc, osc_ref))
    e_dyn_ref = np.array(
        [0.2781795, 0.27817942, 0.28617289, 0.3092472, 0.3092472])
    osc_dyn_ref = np.array([3.87883037e-02, 3.87882995e-02, 7.16099511e-10, 4.46337187e-16,
                            4.46337187e-16])
    e_dyn, osc_dyn = mol.getOscillatorStrengths(dynamic=True)
    assert(np.allclose(e_dyn, e_dyn_ref))
    assert(np.allclose(osc_dyn, osc_dyn_ref))


class ExceptionTests(unittest.TestCase):

    def test_for_data(self):
        mol = Molecule()
        self.assertRaises(Exception, mol.checkData)
        self.assertRaises(Exception, mol.getDFTEnergy)
        self.assertRaises(Exception, mol.getKSTotalEnergy, 1)
        self.assertRaises(Exception, mol.getQPcorrections)
        self.assertRaises(Exception, mol.getOscillatorStrengths)

    def test_existing_molecule_coordinates(self):
        mol = Molecule()
        mol.add_atom("C", 0, 0, 0)
        mol.add_atom("O", 1.3, 0, 0)
        self.assertRaises(Exception, mol.readORB, PATH_TEST / "example.orb")

    def test_existing_molecule_elements(self):
        mol = Molecule()
        mol.add_atom("C", 0, 0, 0)
        mol.add_atom("H", 1.3, 0, 0)
        self.assertRaises(Exception, mol.readORB, PATH_TEST / "example.orb")

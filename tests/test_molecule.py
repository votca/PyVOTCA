from PyVOTCA.molecule import Molecule
from .utils import PATH_TEST
from pathlib import Path


def test_molecule(tmp_path: Path):
    """Check molecule methods."""
    mol = Molecule()

    # check reading method
    mol.readXYZfile(PATH_TEST / "ethylene.xyz")
    assert len(mol.elements) == 6

    # Check writing method
    file_name = tmp_path / "test.xyz"
    mol.writeXYZfile(file_name)
    assert file_name.exists()

"""Test the XTP class."""

import os
from pathlib import Path

from pyvotca import DFTGWBSE, Molecule

from .utils import PATH_TEST


def test_upgrade():
    """Check the mergin between the defauls and the user input."""
    os.environ["VOTCASHARE"] = PATH_TEST.absolute().as_posix()

    # Molecule definition
    mol = Molecule().read_xyz_file(PATH_TEST / "ethylene.xyz")

    user_options = {
        'functional': 'PBE', 'basisset': 'cc-pvtz',
        "dftpackage": {"package": {"name": "orca", "executable": "Path/to/Orca"}},
        "gwbse_engine": {"gwbse_options": {"gwbse": {"mode": 'G0W0'}}}
    }
    file = Path("dftgwbse.xml")
    try:
        votca = DFTGWBSE(mol, options=user_options)
        votca.update_options()
        assert file.exists()
    finally:
        if file.exists():
            os.remove(file)

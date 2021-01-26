"""Test the XTP class."""

import os
from .utils import PATH_TEST
from PyVOTCA import Molecule, XTP


def test_upgrade():
    """Check the mergin between the defauls and the user input."""
    os.environ["VOTCASHARE"] = PATH_TEST.absolute().as_posix()

    # Molecule definition
    mol = Molecule().readXYZfile(PATH_TEST / "ethylene.xyz")

    user_options = {
        'functional': 'PBE', 'basisset': 'cc-pvtz',
        "dftpackage": {"package": {"name": "orca", "executable": "Path/to/Orca"}},
    }
    votca = XTP(mol, options=user_options)

    votca.updateOptions()

    assert False

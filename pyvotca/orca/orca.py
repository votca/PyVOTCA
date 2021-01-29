"""Orca reader/writer module."""
from .molecule import Molecule
from .parsers.orca_parsers import parse_gradient, parse_hessian


class Orca:
    def __init__(self, mol: Molecule):
        self.mol = mol

    def read_gradient(self, gradient_file: str):
        """Read the nuclear gradient from an orca engrad file."""
        # read the gradient from an ORCA engrad calculation
        self.mol.gradient = parse_gradient(gradient_file)

    def read_hessian(self, hessian_file: str):
        """Read Hessian from the orca .hess file."""
        self.mol.hessian = parse_hessian(hessian_file)

"""Orca reader/writer module."""
from .molecule import Molecule
import numpy as np
from .utils import H2EV


class Orca:
    def __init__(self, mol: Molecule):
        self.mol = mol

    def read_gradient(self, gradient_file: str):
        """Reads the nuclear gradient from an orca engrad file."""
        # read the gradient from an ORCA engrad calculation
        natoms = len(self.mol.elements) 
        n = 3 * natoms
        gradient = np.zeros(n)
        count = 0
        nread = 0
        with open(gradfile,'r') as gfile:
            for line in gfile.readlines():
                if '#' not in line:
                    nread+=1
                    if nread > 2:
                        gradient[count] = float(line.split('\n')[0])
                        count+=1
                    if count == n:
                        break

        self.mol.gradient = gradient.reshape(natoms,3)


    def read_hessian(self, hessian_file: str:
        """Reads Hessian from the orca .hess file."""


    def write_gradient(self, gradient_file: str):
        

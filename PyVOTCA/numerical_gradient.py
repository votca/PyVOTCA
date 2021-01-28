"""Numerical gradient."""
import copy
import numpy as np
import os
from .wrapper import XTP
from .molecule import Molecule
import copy as cp

from .utils import BOHR2ANG

__all__ = ["NumericalGradient"]


class NumericalGradient:
    def __init__(self, xtp: XTP, dr: float = 0.001, pathToSimulations='./gradient/'):
        self.xtp = xtp
        self.dr = dr
        self.path = pathToSimulations

    def gen_name(self, name, atom, dir, coord):
        return f"{name}_{atom}_{dir}_{coord}"

    def run_permut(self):
        """ Run's a VOTCA simulation for every displacement of the electric field with strength dE. """
        # Initial structure expected in the mol object of xtp

        # how many atoms
        natoms = len(self.xtp.mol.elements)

        directions = [-1.0, 1.0]
        for atom in range(natoms):
            for coordinate in range(3):
                for direction in directions:
                    # get displaced molecule
                    mol_displaced = Molecule()
                    mol_displaced.copy_and_displace(
                        self.xtp.mol, atom, coordinate, float(direction)*self.dr*BOHR2ANG)
                    name = self.gen_name(
                        mol_displaced.name, atom, direction, coordinate)
                    # make a new xtp wrapper for this one
                    xtp_displaced = XTP(mol_displaced, threads=self.xtp.threads,
                                        options=self.xtp.options, jobname=name, jobdir=self.path)
                    # run this
                    xtp_displaced.run()

    def calcGradient(self, kind, energyLevel=None):
        """ Computes the gradient for a particle/excitation kind expecting all displaced calculation to be available.

        INPUT:  kind of particle/excitation (choices: BSE_singlet, BSE_triplet, QPdiag, QPpert and dft_tot)
                and optionally the energy level if not provided all energy levels will be returned
        OUTPUT: numpy array of nuclear gradient stored in molecule object.     
        """

        # how many atoms
        natoms = len(self.xtp.mol.elements)
        # store gradient in xtp.mol object
        self.xtp.mol.gradient = np.zeros((natoms, 3))

        directions = [-1.0, 1.0]
        for atom in range(natoms):
            for coordinate in range(3):
                E_plus = 0.0
                E_minus = 0.0
                for direction in directions:
                    # get energy for displaced molecules
                    mol_displaced = Molecule()
                    name = self.gen_name(
                        mol_displaced.name, atom, direction, coordinate)
                    orbname = self.path + name + '.orb'
                    mol_displaced.readORB(orbname)
                    if direction > 0:
                        E_plus = mol_displaced.getTotalEnergy(
                            kind, energyLevel)
                    else:
                        E_minus = mol_displaced.getTotalEnergy(
                            kind, energyLevel)

                self.xtp.mol.gradient[atom, coordinate] = (
                    E_plus - E_minus)/(2.0 * self.dr)

        self.xtp.mol.hasGradient = True
        return self.xtp.mol.gradient

    def getGradient(self, kind, energyLevel=None):
        self.run_permut()
        self.calcGradient(kind, energyLevel)

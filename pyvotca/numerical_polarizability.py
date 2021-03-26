"""Numerical polarizability.

API
---
.. autoclass:: NumericalPolarizability

"""

import numpy as np
import os
import itertools
import h5py

from .molecule import Molecule
from .utils import BOHR2ANG
from .xtp import DFTGWBSE

__all__ = ["NumericalPolarizability"]


class NumericalPolarizability:
    """Compute the polarizability numerically using XTP."""
    def __init__(self,
                 xtp: DFTGWBSE,
                 dE: float = 0.002,
                 path_to_simulations: str = './polar_experiments/'):
        self.xtp = xtp
        self.dE = dE
        self.path = path_to_simulations

    def mkfldr(self, path: str):
        try:
            os.mkdir(path)
        except OSError:
            print("Creation of the directory %s failed" % path)

    def gen_name(self, name, E) -> str:
        return f"{name}_E{E[0]:.4f}E{E[1]:.4f}E{E[2]:.4f}"

    def run_permut(self):
        """Run a VOTCA simulation for every displacement of the electric field with displacement dE in atomic units."""
        # compute all necessary directions of the field (permutations)
        b = [self.dE, -self.dE, 0.0]
        a = [b, b, b]
        perm = [element for element in itertools.product(*a)]
        perm = [s for s in perm if all(s) == False]
        # create a folder to contain all the results from the experiments
        if (not os.path.exists(self.path)):
            self.mkfldr(self.path)
        counter = 1
        for E in perm:
            self.xtp.options.dftpackage = {
                "package": {
                    "use_external_field": "true",
                    "externalfield": f"{E[0]} {E[1]} {E[2]}"
                }
            }
            name = self.gen_name(self.xtp.mol.name, E)
            xtp_withField = DFTGWBSE(self.xtp.mol,
                                     threads=self.xtp.threads,
                                     options=self.xtp.options,
                                     jobname=name,
                                     jobdir=self.path)
            print(
                f"Running XTP for field: X{E[0]:7.3f} Y{E[1]:7.3f} Z{E[2]:7.3f}  {counter}/{len(perm)}"
            )
            xtp_withField.run()
            counter += 1

    def getEnergiesFromOrbFile(self, kind: str, i: int, j: int, k: int):
        """ Read the energies from the orb file for a certain kind of particle/excitation.

        i,j,k represents the direction of the electric field 1,0,0 corresponds to dE,0,0 etc.
        OUTPUT: numpy array with the energies of the particle/excitation kind.
        """
        nameOrbFile = self.path + self.gen_name(
            self.xtp.mol.name, np.array([i * self.dE, j * self.dE, k * self.dE
                                         ])) + ".orb"
        f = h5py.File(nameOrbFile, 'r')
        orb = f['QMdata']
        if (kind == 'BSE_singlet' or kind == 'BSE_triplet' or kind == 'QPdiag'):
            return (np.array(orb[kind]['eigenvalues'][()]))
        elif (kind == 'QPpert'):
            return (np.array(orb['QPpert_energies'][()]))
        elif (kind == 'dft_tot'):
            return (np.array(orb.attrs['qm_energy']))
        else:
            raise Exception(
                'Invalid kind, kind of particle/excitation should be one of:' +
                ' BSE_singlet, BSE_triplet, QPdiag, QPpert and dft_tot\n' +
                f'the value of kind was: {kind}')

    def calc_polarizability(self, kind: str, energyLevel=None) -> np.ndarray:
        """ Computes the polarizability tensor for a particle/excitation type.

        INPUT:  kind of particle/excitation (choices: BSE_singlet, BSE_triplet,
                QPdiag, QPpert and dft_tot) and optionally the energy level if
                not provided all energy levels will be returned 
        OUTPUT: numpy array of polarizability tensors.
        """
        # Get all the energies from the files (encoding mop = -1 0 1)
        ooo = self.getEnergiesFromOrbFile(kind, 0, 0, 0)
        moo = self.getEnergiesFromOrbFile(kind, -1, 0, 0)
        poo = self.getEnergiesFromOrbFile(kind, 1, 0, 0)
        omo = self.getEnergiesFromOrbFile(kind, 0, -1, 0)
        opo = self.getEnergiesFromOrbFile(kind, 0, 1, 0)
        oom = self.getEnergiesFromOrbFile(kind, 0, 0, -1)
        oop = self.getEnergiesFromOrbFile(kind, 0, 0, 1)

        mmo = self.getEnergiesFromOrbFile(kind, -1, -1, 0)
        pmo = self.getEnergiesFromOrbFile(kind, 1, -1, 0)
        mpo = self.getEnergiesFromOrbFile(kind, -1, 1, 0)
        ppo = self.getEnergiesFromOrbFile(kind, 1, 1, 0)

        mom = self.getEnergiesFromOrbFile(kind, -1, 0, -1)
        pom = self.getEnergiesFromOrbFile(kind, 1, 0, -1)
        mop = self.getEnergiesFromOrbFile(kind, -1, 0, 1)
        pop = self.getEnergiesFromOrbFile(kind, 1, 0, 1)

        omm = self.getEnergiesFromOrbFile(kind, 0, -1, -1)
        opm = self.getEnergiesFromOrbFile(kind, 0, 1, -1)
        omp = self.getEnergiesFromOrbFile(kind, 0, -1, 1)
        opp = self.getEnergiesFromOrbFile(kind, 0, 1, 1)

        # Compute second derivatives of the energy with respect to the field
        # strength
        coeff = (-1. / (self.dE**2))
        axx = coeff * (moo - 2 * ooo + poo)
        ayy = coeff * (omo - 2 * ooo + opo)
        azz = coeff * (oom - 2 * ooo + oop)
        axy = 0.25 * coeff * (mmo - pmo - mpo + ppo)
        axz = 0.25 * coeff * (mom - pom - mop + pop)
        ayz = 0.25 * coeff * (omm - opm - omp + opp)

        # Create polarization tensor
        polar_m = np.squeeze(np.array([[axx, axy, axz], [axy, ayy, ayz],
                                       [axz, ayz, azz]]),
                             axis=-1)

        # if the entries of the matrix are arrays we need to reshape to get an
        # array of matrices
        if (axx.size > 1):
            polar_m = np.moveaxis(polar_m, [2], [0])
        else:  # there can't be an energyLevel
            energyLevel = None

        if (energyLevel is None):
            return polar_m
        else:
            return polar_m[energyLevel]

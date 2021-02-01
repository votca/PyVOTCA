"""Electron-phonon-coupling module."""
from .molecule import Molecule
import numpy as np
import periodictable
import matplotlib.pyplot as plt
from .utils import H2EV
from .utils import AFU2EV
import scipy.linalg as la
from typing import List, Optional, Tuple, Union


class Electronphonon:
    def __init__(self): pass

    def get_mass(self, elements):
        mass = []
        for atom in elements:
            mass.append(periodictable.formula(atom).mass)
        mass = np.array(mass, dtype=float)
        return mass

    def calculate_electron_phonon_couplings(self, elements, gs_hessian, gradient, plot=False) -> Tuple[np.ndarray, np.ndarray]:
        """Calculates the mode-resolved electron-phonon coupling from 
           the normal mode projected excited state gradient:

                V_nu^ep = 1/(M * w_nu^2) * | e_v * gradient(E) |^2

            where:
                - w_nu: frequency of the nu-th vibrational mode in GS
                - e_nu: eigenvector of the nu-th vibrational mode in GS

            INPUT: - elements:  np.array with N chemical elements 
                   - hessian:   np.array of GS Hessian (3N,3N)
                                in (a.f.u)^2 = Hartree/Bohr^2/amu
                   - gradient:  flattened excited state gradient (3N)
                                in Hartree/Bohr

            OUTPUT: - freq:         vibrational frequencies in a.f.u
                    - ep_couplings: electron-phonon coupling per mode in Hartree
        """

        # we don't want to change the input gradient
        es_gradient = np.copy(gradient)

        # get vector with masses of the atoms
        mass = np.repeat(self.get_mass(elements), 3)

        # get vibrational frequencies and eigenmodes
        freq, modes = self.calculate_vibrational_modes(elements, gs_hessian)

        # flatten the excited state gradient, if necessary
        es_gradient = es_gradient.flatten()

        # get the normal-mode projected gradient (already divided by sqrt(M))
        es_gradient /= np.sqrt(mass)
        nm_gradient = modes.T.dot(es_gradient)

        # ignore any contribution from translation/rotation and imaginary frequencies
        nm_gradient[np.where(freq < 1e-5)] = 0.

        # electron-phonon couplings in Hartree
        ep_couplings = (nm_gradient/freq)**2

        if plot:
            plt.bar(AFU2EV*freq, H2EV*ep_couplings, width=0.005, alpha=0.6)
            plt.plot(AFU2EV*freq, H2EV*np.cumsum(ep_couplings), lw=2., alpha=0.8)
            plt.ylabel('cum. electron-phonon coupling (eV)')
            plt.xlabel('hw (eV)')
            plt.show()

        return (freq, ep_couplings)

    def calculate_vibrational_modes(self, elements, hessian_in):
        """Calculates the vibrational frequencies and eigenmodes from a Hessian matrix and a list of elements"""

        # we don't want to change the input hessian
        hessian = np.copy(hessian_in)

        # get vector with masses of the atoms
        mass = np.repeat(self.get_mass(elements), 3)

        # mass weighting the Hessian
        hessian /= np.sqrt(mass)
        hessian /= np.sqrt(mass.reshape((len(mass), 1)))

        # get eigenmodes (squared frequency in (a.f.u)^2 = Hartree/Bohr^2/amu)
        freq_sq, modes = la.eigh(hessian)

        # convert to frequencies in a.f.u
        cfreq = np.array(freq_sq, dtype=np.complex128)
        cfreq = np.sqrt(cfreq)

        # check for imaginary frequencies and store them as negative real ones
        freq = np.where(np.isreal, np.real(cfreq), -np.imag(cfreq))

        # sort frequencies and eigenvectors
        isort = freq.argsort()
        modes = np.real(modes[:, isort])
        freq = np.sort(freq)

        return (freq, modes)

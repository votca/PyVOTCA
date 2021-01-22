"""Visualization module."""
from .wrapper import XTP
import numpy as np
import matplotlib.pyplot as plt
from .utils import Utils


class Visualization:
    def __init__(self, votca: XTP):
        self.votca = votca
        self.constants = Utils()

    def plotQPcorrections(self):
        QPcorrections = self.votca.getQPcorrections()
        qpmin = self.votca.qpmin
        qpmax = self.votca.qpmax + 1
        correctedKS = self.votca.KSenergies[qpmin:qpmax]
        plt.plot(correctedKS, QPcorrections, 'ro')
        plt.xlabel('KS Energy (eV)')
        plt.ylabel('QP correction (eV)')
        plt.show()

    def plotAbsorptionGaussian(self, dynamic=False, min=0.0, max=10.0, points=1000, sigma=0.2):

        energy, osc = self.votca.getOscillatorStrengths(dynamic)

        # convert energies from Hartree to eV
        energy *= self.constants.h2ev  # to eV
        # make a stick plot with oscillator strength
        plt.stem(energy, osc, basefmt=" ")
        # apply unormalized Gaussian lineshape
        e = np.linspace(min, max, points)
        spectrum = 0
        for i in range(len(energy)):
            spectrum += osc[i] * self.gaussian(e, energy[i], sigma)

        plt.plot(e, spectrum, 'k', linewidth=2)
        plt.ylim(bottom=0)
        plt.title(f'Gaussian lineshape with sigma = {sigma}eV')
        plt.xlabel('Energy (eV)')
        plt.ylabel('Absorption (arb. units)')
        plt.show()

    def gaussian(self, x, mu, sig):
        #ATTN: not normalized
        return np.exp(-0.5 * ((x - mu) / sig) ** 2)

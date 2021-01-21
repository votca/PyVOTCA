"""Visualization module."""
from .wrapper import XTP
import numpy as np
import matplotlib.pyplot as plt


class Visualization:
    def __init__(self, votca: XTP):
        self.votca = votca

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
        # get energies/oscillator strengths
        if dynamic:
            energy = self.votca.BSE_singlet_energies_dynamic
        else:
            energy = self.votca.BSE_singlet_energies
        td = self.votca.transition_dipoles
        osc = []
        for i in range(len(energy)):
            osc.append(2./3. * energy[i] * np.sum(np.power(td[i],2)))

        # convert energies from Hartree to eV
        energy *= 27.21138 # to eV
        # make a stick plot with oscillator strength
        plt.stem(energy,osc,basefmt=" ")
        # apply unormalized Gaussian lineshape
        e = np.linspace(min,max,points)
        spectrum = 0
        for i in range(len(energy)):
            spectrum += osc[i] * self.gaussian(e,energy[i], sigma)

        plt.plot(e,spectrum,'k',linewidth=2)
        plt.ylim(bottom=0)
        plt.title(f'Gaussian lineshape with sigma = {sigma}eV')
        plt.xlabel('Energy (eV)')
        plt.ylabel('Absorption (arb. units)')
        plt.show()

    def gaussian(self, x, mu, sig):
        #ATTN: not normalized
        return np.exp(-np.power((x - mu)/sig, 2.)/2)
        

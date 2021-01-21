import numpy as np
import os
from .wrapper import XTP
import matplotlib.pyplot as plt

class Visualization:
    def __init__(self, votca : XTP ): 
        self.votca = votca

    def plotQPcorrections(self):
        QPcorrections = self.votca.getQPcorrections()
        qpmin = self.votca.qpmin
        qpmax = self.votca.qpmax+1
        correctedKS = self.votca.KSenergies[qpmin:qpmax]
        plt.plot(correctedKS,QPcorrections,'ro')
        plt.xlabel('KS Energy (eV)')
        plt.ylabel('QP correction (eV)')
        plt.show()



"""Numerical gradient."""
import copy
import numpy as np
import os
from collections import namedtuple
from itertools import cycle


#from .wrapper import VOTCA as votca

__all__ = ["NumericalGradient"]

"""**************************************************************
* PART I: Functions to run a simulation for all electric fields *
**************************************************************"""


def mkfldr(path):
    try:
        os.mkdir(path)
    except OSError:
        print("Creation of the directory %s failed" % path)


"""**********************************************************
* PART II: Gradient Calculator                              *
**********************************************************"""


class NumericalGradient:
    def __init__(self, dr, n_atoms, name="",  pathToSimulations='./experiments'):
        self.dr = dr
        self.n_atoms = n_atoms
        self.name = name
        self.path = pathToSimulations

    def run_permut(self, dr, name, threads: int = 4):
        """ Run's a VOTCA simulation for every direction (perm) of the electric field with strength dE. """
        # Read XYZ
        filename = "{}.xyz".format(str(name))
        xyzfile = open(filename)
        xyz = cio.read_xyz(xyzfile)
        xyzfile.close()

        # create a folder to contain all the results from the different experiments
        if(not os.path.exists('./experiments')):
            mkfldr('./experiments')
        # run VOTCA for every displacement
        counter = 1
        for atom in range(int(xyz.coords.size / 3)):
            for coordinate in range(3):
                # displacement plus
                xyz_plus = copy.deepcopy(xyz)
                xyz_plus.coords[atom, coordinate] += dr
                xyzfilename_plus = "{}_at{}_dir{}_plus.xyz".format(
                    str(name), str(atom), str(coordinate))
                xyzfile_plus = open(xyzfilename_plus, "w")
                cio.write_xyz(xyzfile_plus, *xyz_plus)
                xyzfile_plus.close()
                print("Running XTP for displacement +{} for atom {} in direction {}".format(
                    str(dr), str(atom), str(coordinate)))
                votca.run(atom, coordinate, name, "plus", threads)

                # displacement minus
                xyz_minus = copy.deepcopy(xyz)
                xyz_minus.coords[atom, coordinate] -= dr
                xyzfilename_minus = "{}_at{}_dir{}_minus.xyz".format(
                    str(name), str(atom), str(coordinate))
                xyzfile_minus = open(xyzfilename_minus, "w")
                cio.write_xyz(xyzfile_minus, *xyz_minus)
                xyzfile_minus.close()
                print("Running XTP for displacement -{} for atom {} in direction {}".format(
                    str(dr), str(atom), str(coordinate)))
                votca.run(atom, coordinate, name, "minus", threads)

    

    def getGradient(self, kind, energyLevel=None):
        """ Computes the gradient for a particle/excitation kind.

        INPUT:  kind of particle/excitation (choices: BSE_singlet, BSE_triplet, QPdiag, QPpert and dft_tot)
                and optionally the energy level if not provided all energy levels will be returned
        OUTPUT: numpy array of polarizability tensors.     
        """

        print(self.n_atoms)

        gradient = np.zeros((self.n_atoms, 3))

        for atom in range(self.n_atoms):
            for coordinate in range(3):
                # get energies from files
                filename = "./experiments/{}_at{}_dir{}_{}.orb".format(
                    str(self.name), str(atom), str(coordinate), "plus")
                Eplus = votca.getEnergies(
                    kind, filename, energyLevel)
                filename = "./experiments/{}_at{}_dir{}_{}.orb".format(
                    str(self.name), str(atom), str(coordinate), "minus")
                Eminus = votca.getEnergies(
                    kind, filename, energyLevel)

                # Compute derivative
                gradient[atom, coordinate] = (
                    Eplus - Eminus)/(2.0*self.dr*1.8897259886)

        return gradient

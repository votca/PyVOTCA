"""Molecule representation."""
from pathlib import Path
from typing import Optional, Union

import h5py
import numpy as np

from .utils import BOHR2ANG

Pathlike = Union[Path, str]


class Molecule:
    """Molecule definition."""

    def __init__(self):
        self.elements = []
        self.coordinates = []
        self.name = "molecule"
        self.hasData = False
        self.hasXYZ = False

    def add_atom(self, element: str, x: float, y: float, z: float):
        self.elements.append(element)
        self.coordinates.append(np.array([x, y, z]))
        self.hasXYZ = True

    def printXYZ(self):
        for (element, coordinates) in zip(self.elements, self.coordinates):
            print(element, coordinates)

    def readXYZfile(self, filename: Pathlike):
        with open(filename, 'r') as handler:
            lines = handler.readlines()

        self.name = Path(filename).stem
        arr = [(row[0], np.array(row[1:], dtype=float)) for row in [
            x.split() for x in lines[2:]]]
        self.elements, self.coordinates = tuple(zip(*arr))
        self.hasXYZ = True

    def writeXYZfile(self, filename: Pathlike):
        """Write the molecule in XYZ format."""
        atoms = "\n".join(f"{elem} {xyz[0]:.4f} {xyz[1]:.4f} {xyz[2]:.4f}" for elem, xyz in zip(
            self.elements, self.coordinates))
        mol = f"""{len(self.elements)}
{self.name} created by PyVOTCA writer
{atoms}
"""
        with open(filename, "w") as xyzfile:
            xyzfile.write(mol)

    def getDFTEnergy(self):
        """ Returns the DFT total energy."""
        self.checkData()

        return self.DFTenergy

    def getKSTotalEnergy(self, level=''):
        """Returns the excited state KS total energy."""
        self.checkData()

        lumo = self.homo + 1

        total_energy = self.DFTenergy
        if (level < lumo):
            return(total_energy - self.KSenergies[level])
        elif level < len(self.KSenergies):
            return(total_energy + self.KSenergies[level])
        else:
            print("Requested KS level {} does not exist.")
            return 0.0

    def getQPTotalEnergy(self, level=''):
        """Returns the excited state QP total energy."""
        self.checkData()

        lumo = self.homo + 1

        total_energy = self.DFTenergy
        if (level < lumo):
            return(total_energy - self.QPenergies[level - self.qpmin])
        elif level < len(self.KSenergies):
            return(total_energy + self.QPenergies[level - self.qpmin])
        else:
            print("Requested QP level {} does not exist.")
            return 0.0

    def getQPdiagTotalEnergy(self, level=''):
        """Returns the excited state diag QP total energy."""
        self.checkData()

        lumo = self.homo + 1

        total_energy = self.DFTenergy
        if (level < lumo):
            return(total_energy - self.QPenergies_diag[level - self.qpmin])
        elif level < len(self.KSenergies):
            return(total_energy + self.QPenergies_diag[level - self.qpmin])
        else:
            print("Requested diag QP level {} does not exist.")
            return 0.0

    def getBSEsingletTotalEnergy(self, level=''):
        """ Returns the excited state BSE Singlet total energy."""
        self.checkData()

        if level < len(self.BSE_singlet_energies):
            return(self.DFTenergy + self.BSE_singlet_energies[level])
        else:
            print("Requested BSE singlet level {} does not exist.")
            return 0.0

    def getBSEtripletTotalEnergy(self, level=''):
        """ Returns the excited state BSE Singlet total energy."""
        self.checkData()

        if level < len(self.BSE_triplet_energies):
            return(self.DFTenergy + self.BSE_triplet_energies[level])
        else:
            print("Requested BSE triplet level {} does not exist.")
            return 0.0

    def getBSEsingletDynamicTotalEnergy(self, level=''):
        """ Returns the excited state BSE Singlet total energy."""
        self.checkData()

        if level < len(self.BSE_singlet_energies_dynamic):
            return(self.DFTenergy + self.BSE_singlet_energies_dynamic[level])
        else:
            print("Requested dynamic BSE singlet level {} does not exist.")
            return 0.0

    def getBSEtripletDynamicTotalEnergy(self, level=''):
        """ Returns the excited state BSE Singlet total energy."""
        self.checkData()

        if level < len(self.BSE_triplet_energies_dynamic):
            return(self.DFTenergy + self.BSE_triplet_energies_dynamic[level])
        else:
            print("Requested dynamic BSE triplet level {} does not exist.")
            return 0.0

    # Parse energies/info from HDF5

    def readORB(self, orbfile):

        with h5py.File(orbfile, 'r') as handler:
            orb = handler['QMdata']
            # get coordinates
            atoms = orb['qmmolecule']['qmatoms']
            # coordinates are stored in Bohr!
            arr = [(atom['element'][0].decode(), BOHR2ANG * np.array(
                [atom['posX'][0], atom['posY'][0], atom['posZ'][0]], dtype=float)) for atom in atoms]
            elements_in, coordinates_in = tuple(zip(*arr))

            if not self.hasXYZ:
                self.elements = elements_in
                self.coordinates = coordinates_in
            else:
                for i in range(len(elements_in)):
                    if not self.elements[i] == elements_in[i]:
                        raise Exception('Element {} in molecule differs from element {} in orb file!'.format(
                            self.elements[i], elements_in[i]))

                    if not np.allclose(self.coordinates[i], coordinates_in[i]):
                        raise Exception('Coordinates of element {} in molecle {} differ from coordinates in orb file {}'.format(
                            i, self.coordinates[i], coordinates_in[i]))

            self.hasXYZ = True

            self.homo = int(orb.attrs['occupied_levels']) - 1
            self.DFTenergy = float(orb.attrs['qm_energy'])
            self.KSenergies = read_flatten_array(orb, 'mos', 'eigenvalues')
            self.QPenergies = read_flatten_array(orb, 'QPpert_energies')
            self.QPenergies_diag = read_flatten_array(
                orb, 'QPdiag', 'eigenvalues')
            self.BSE_singlet_energies = read_flatten_array(
                orb, 'BSE_singlet', 'eigenvalues')
            self.BSE_triplet_energies = read_flatten_array(
                orb, 'BSE_triplet', 'eigenvalues')
            self.BSE_singlet_energies_dynamic = read_flatten_array(
                orb, 'BSE_singlet_dynamic')
            self.BSE_triplet_energies_dynamic = read_flatten_array(
                orb, 'BSE_triplet_dynamic')
            self.qpmin = int(orb.attrs['qpmin'])
            self.qpmax = int(orb.attrs['qpmax'])
            td = orb['transition_dipoles']
            self.transition_dipoles = np.array(
                [td[dset][()] for dset in td.keys()])
            self.hasData = True

    def getQPcorrections(self):
        self.checkData()

        QPcorrections = self.QPenergies - \
            self.KSenergies[self.qpmin:self.qpmin + len(self.QPenergies)]

        return QPcorrections.flatten()

    def getOscillatorStrengths(self, dynamic=False):
        self.checkData()

        # get energies/oscillator strengths
        if dynamic:
            energy = self.BSE_singlet_energies_dynamic
        else:
            energy = self.BSE_singlet_energies
        osc = [(2. / 3.) * e * (t ** 2).sum()
               for e, t in zip(energy, self.transition_dipoles)]

        return energy, np.array(osc)

    def checkData(self):
        if not self.hasData:
            raise Exception("No energy has been stored!")


def read_flatten_array(group: h5py.Group, key1: str, key2: Optional[str] = None):
    """Read an array from h5py handler and flatten it."""
    if key2 is None:
        arr = group[key1][()]
    else:
        arr = group[key1][key2][()]

    return arr.flatten()

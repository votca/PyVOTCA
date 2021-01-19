import numpy as np
import os
import h5py
import xml.etree.ElementTree as ET
import itertools
import copy
from itertools import cycle
from collections import namedtuple

""" XYZ FILES """
def read_xyz(fin):
    """ read a xyz file from file handle
    Parameters
    ----------
    fin : file handle
        file to read from
    Returns
    -------
    fin : open file
    xyz : namedtuple
        returns a named tuple with coords, title and list of atomtypes.
    See Also
    --------
    write_xyz
    """
    natoms = int(fin.readline())
    title = fin.readline()[:-1]
    coords = np.zeros([natoms, 3], dtype="float64")
    atomtypes = []
    for x in coords:
        line = fin.readline().split()
        atomtypes.append(line[0])
        x[:] = list(map(float, line[1:4]))

    return namedtuple("XYZFile", ["coords", "title", "atomtypes"]) \
        (coords, title, atomtypes)


def write_xyz(fout, coords, title="", atomtypes=("A",)):
    """ write a xyz file from file handle
    Writes coordinates in xyz format. It uses atomtypes as names. The list is
    cycled if it contains less entries than there are coordinates,
    One can also directly write xyz data which was generated with read_xyz.
    >>> xx = read_xyz("in.xyz")
    >>> write_xyz(open("out.xyz", "w"), *xx)
    Parameters
    ----------
    fout : an open file
    coords : np.array
        array of coordinates
    title : title section, optional
        title for xyz file
    atomtypes : iteratable
        list of atomtypes.
    See Also
    --------
    read_xyz
    """
    fout.write("%d\n%s\n" % (coords.size / 3, title))
    for x, atomtype in zip(coords.reshape(-1, 3), cycle(atomtypes)):
        fout.write("%s %.18g %.18g %.18g\n" % (atomtype, x[0], x[1], x[2]))


"""**************************************************************
* PART I: Functions to run a simulation for all electric fields *
**************************************************************"""


def mkfldr(path):
    try:
        os.mkdir(path)
    except OSError:
        print("Creation of the directory %s failed" % path)


def run_votca(atom, direction, name, case, threads):
    """ Runs VOTCA and moves results to the experiments folder """
    name_disp = "{}_at{}_dir{}_{}".format(str(name), str(atom), str(direction),str(case))
    votcacmd = "xtp_tools -e dftgwbse -o dftgwbse.xml -n {} -t {} > ./experiments/logfile_{}.log".format(str(name_disp), str(threads), str(name_disp))
    os.system(votcacmd)
    os.replace('{}.orb'.format(str(name_disp)), './experiments/{}.orb'.format( str(name_disp)))





            

"""**********************************************************
* PART II: Gradient Calculator                              *
**********************************************************"""
class NumericalGradient:
    def __init__(self, dr, n_atoms, name,  pathToSimulations = './experiments'):
        self.dr = dr
        self.n_atoms = n_atoms
        self.name = name
        self.path = pathToSimulations 


    def run_permut(dr, name, threads):
        """ Run's a VOTCA simulation for every direction (perm) of the electric field with strength dE. """
        # Read XYZ
        filename = "{}.xyz".format(str(name))
        xyzfile = open(filename)
        xyz=read_xyz(xyzfile)
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
                xyz_plus.coords[atom,coordinate]+=dr
                xyzfilename_plus = "{}_at{}_dir{}_plus.xyz".format(str(name),str(atom),str(coordinate))
                xyzfile_plus=open(xyzfilename_plus, "w")
                write_xyz(xyzfile_plus,*xyz_plus)
                xyzfile_plus.close()
                print("Running XTP for displacement +{} for atom {} in direction {}".format(str(dr),str(atom),str(coordinate)))
                run_votca(atom, coordinate, name, "plus", threads)

                # displacement minus
                xyz_minus = copy.deepcopy(xyz)
                xyz_minus.coords[atom,coordinate]-=dr
                xyzfilename_minus = "{}_at{}_dir{}_minus.xyz".format(str(name),str(atom),str(coordinate))
                xyzfile_minus=open(xyzfilename_minus, "w")
                write_xyz(xyzfile_minus,*xyz_minus)
                xyzfile_minus.close()
                print("Running XTP for displacement -{} for atom {} in direction {}".format(str(dr),str(atom),str(coordinate)))
                run_votca(atom, coordinate, name, "minus", threads)





        
    def getEnergiesFromOrbFile(self, kind, filename, level):
        """ Read the energies from the orb file for a certain kind of particle/excitation.

        OUTPUT: numpy array with the energies of the particle/excitation kind.
        """
        nameOrbFile = filename
        f = h5py.File(nameOrbFile, 'r')
        orb = f['QMdata']
        occupied_levels = int(orb.attrs['occupied_levels'])
        
        total_energy = orb.attrs['qm_energy']
        if(kind == 'BSE_singlet' or kind == 'BSE_triplet'):
            return(total_energy+orb[kind]['eigenvalues'][level])
        elif (kind == 'QPdiag'):
            if (level < occupied_levels):
                return(total_energy - orb[kind]['eigenvalues'][level])
            else:
                return(total_energy + orb[kind]['eigenvalues'][level])
        elif (kind == 'QPpert'):
            if (level < occupied_levels):
                return(total_energy - orb['QPpert_energies'][level])
            else:
                return(total_energy + orb['QPpert_energies'][level])
        elif (kind == 'dft_tot'):
            return total_energy
        else:
            print("Invalid kind!")
            exit()

    def getGradient(self, kind, energyLevel=None):
        """ Computes the gradient for a particle/excitation kind.

        INPUT:  kind of particle/excitation (choices: BSE_singlet, BSE_triplet, QPdiag, QPpert and dft_tot)
                and optionally the energy level if not provided all energy levels will be returned
        OUTPUT: numpy array of polarizability tensors.     
        """

        print(self.n_atoms)

        gradient=np.zeros((self.n_atoms,3))

        for atom in range(self.n_atoms):
            for coordinate in range(3):
                # get energies from files
                filename = "./experiments/{}_at{}_dir{}_{}.orb".format(str(self.name), str(atom), str(coordinate),"plus")
                Eplus=self.getEnergiesFromOrbFile(kind,filename, energyLevel)
                filename = "./experiments/{}_at{}_dir{}_{}.orb".format(str(self.name), str(atom), str(coordinate),"minus")
                Eminus=self.getEnergiesFromOrbFile(kind,filename, energyLevel)

                # Compute derivative
                gradient[atom,coordinate] = (Eplus - Eminus)/(2.0*self.dr*1.8897259886)


        return gradient
 
    

import os
import h5py
import numpy as np
import xml.etree.ElementTree as ET
from .molecule import Molecule


__all__ = ["XTP"]


class XTP:

    def __init__(self, threads=1, jobname='dftgwbse', jobdir='./'): 
        self.threads = threads
        self.jobname = jobname
        self.jobdir = jobdir
        self.orbfile = ''
        self.hasData = False
        self.options = dict() 


    def updateOptions(self):

        # parsing defaults 
        votcashare=os.environ.get('VOTCASHARE')
        default_options=votcashare+'/xtp/xml/dftgwbse.xml'
        options = ET.parse(default_options) 
        root = options.getroot()

        for option, value in self.options.items():
            if value != '':
                for setting in root.iter(option):
                    setting.text = str(value)

        # write out xml
        options.write('dftgwbse.xml') 


    # just runs xtp_tools with CMDline call
    def run(self, mol : Molecule):
        # update and write the options
        self.updateOptions()

        # write the XYZfile
        xyzname=mol.getName()
        xyzfile=xyzname+".xyz"
        mol.writeXYZfile(xyzfile)

        """ Runs VOTCA and moves results a job folder, if requested """
        if not os.path.exists(self.jobdir):
            os.makedirs(self.jobdir)

        votcacmd = "xtp_tools -e dftgwbse -o dftgwbse.xml -n {} -t {} > {}{}.log".format(str(xyzname), str(self.threads), str(self.jobdir), str(self.jobname))
        os.system(votcacmd)
        # copy orbfile, if jobdir is not default
        if (self.jobdir != "./"):
            self.orbfile = '{}{}.orb'.format(str(self.jobdir),str(xyzname))
            os.replace('{}.orb'.format(str(xyzname)),self.orbfile)
        else:
            self.orbfile = '{}.orb'.format(str(xyzname))
        
        self.getEnergies(self.orbfile)


    # Reads energies from an existing HDF5 orb file
    def getTotalEnergy(self, kind, level, orbfile=''):
        """ Read the energies from the orb file for a certain kind of particle/excitation.

        OUTPUT: numpy array with the energies of the particle/excitation kind.
        """
        if self.hasData == False:
            print("No energy has been stored!")
            exit(0)

        occupied_levels = self.homo

        total_energy = self.DFTenergy
        if (kind == 'BSE_singlet'):
           return(total_energy + self.BSE_singlet_energies[level])
        elif (kind == 'BSE_triplet'):
            return(total_energy + self.BSE_triplet_energies[level])
        elif (kind == 'QPdiag'):
            if (level < occupied_levels):
                return(total_energy - self.QPenergies_diag[level-qpmin])
            else:
                return(total_energy + self.QPenergies_diag[level-qpmin])
        elif (kind == 'QPpert'):
            if (level < occupied_levels):
                return(total_energy - self.QPenergies[level-qpmin])
            else:
                return(total_energy + self.QPenergies[level-qpmin])
        elif (kind == 'dft_tot'):
            return total_energy
        else:
            print("Invalid kind!")
            exit()


    # Parse energies/info from HDF5
    def getEnergies(self,orbfile=''):
        nameOrbFile=''
        if (self.orbfile == '' ):
            if (orbfile == ''):
                print("No HDF5 file for reading is known")
                exit(0)
            else:
                nameOrbFile = orbfile
        else:
            nameOrbFile = self.orbfile

        f = h5py.File(nameOrbFile, 'r')
        orb = f['QMdata']
        self.homo = int(orb.attrs['occupied_levels'])
        self.DFTenergy = float(orb.attrs['qm_energy'])
        self.KSenergies = np.array(orb['mos']['eigenvalues'][:])
        self.QPenergies = np.array(orb['QPpert_energies'][:])
        self.QPenergies_diag = np.array(orb['QPdiag']['eigenvalues'][:])
        self.BSE_singlet_energies = np.array(orb['BSE_singlet']['eigenvalues'][:])
        self.BSE_triplet_energies = np.array(orb['BSE_triplet']['eigenvalues'][:])
        self.BSE_singlet_energies_dynamic = np.array(orb['BSE_singlet_dynamic'][:])
        self.BSE_triplet_energies_dynamic = np.array(orb['BSE_triplet_dynamic'][:])
        self.qpmin = int(orb.attrs['qpmin'])
        self.qpmax = int(orb.attrs['qpmax'])
        self.hasData = True
        

    def getKSenergies(self):
        return self.KSenergies


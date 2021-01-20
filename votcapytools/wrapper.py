import os
import h5py

__all__ = ["VotcaWrapper"]


class VotcaWrapper:

    def __init__(self): pass

    # just runs xtp_tools with CMDline call, expects
    # - an xyz file xyzname.xyz
    # - the dftgwbse.xml options file
    # is present in the same directory
    def run(self, xyzname, threads=1, jobname="dftgwbse", jobdir="./"):
        """ Runs VOTCA and moves results a job folder, if requested """
        #name_disp = "{}_at{}_dir{}_{}".format(str(name), str(atom), str(direction), str(case))
        votcacmd = "xtp_tools -e dftgwbse -o dftgwbse.xml -n {} -t {} > {}{}.log".format(str(xyzname), str(threads), str(jobdir), str(jobname))
        os.system(votcacmd)
        # copy orbfile, if jobdir is not default
        if (jobdir != "./"):
            os.replace('{}.orb'.format(str(xyzname)),'{}{}.orb'.format(str(jobdir),str(xyzname)))

    # Reads energies from an existing HDF5 orb file
    def getEnergies(self, kind, filename, level):
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

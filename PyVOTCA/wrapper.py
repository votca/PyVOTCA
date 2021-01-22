"""Votca Wrapper."""
import os
import subprocess
import xml.etree.ElementTree as ET


import numpy as np

from .molecule import Molecule

__all__ = ["XTP"]


class XTP:

    def __init__(self, mol: Molecule, threads=1, jobname='dftgwbse', jobdir='./'):
        self.mol = mol
        self.threads = threads
        self.jobname = jobname
        self.jobdir = jobdir
        self.orbfile = ''
        self.options = {}

    def updateOptions(self):

        # parsing defaults
        votcashare = os.environ.get('VOTCASHARE')
        default_options = f'{votcashare}/xtp/xml/dftgwbse.xml'
        options = ET.parse(default_options)
        root = options.getroot()

        for option, value in self.options.items():
            if value != '':
                for setting in root.iter(option):
                    setting.text = str(value)

        # write out xml
        options.write('dftgwbse.xml')

    # just runs xtp_tools with CMDline call

    def run(self):
        # update and write the options
        self.updateOptions()

        # write the XYZfile
        xyzname = self.mol.name
        xyzfile = xyzname + ".xyz"
        self.mol.writeXYZfile(xyzfile)

        """ Runs VOTCA and moves results a job folder, if requested """
        if not os.path.exists(self.jobdir):
            os.makedirs(self.jobdir)

        votcacmd = f"xtp_tools -e dftgwbse -o dftgwbse.xml -n {xyzname} -t {self.threads} > {self.jobdir}{self.jobname}.log"

        subprocess.check_output(votcacmd, shell=True, stderr=subprocess.STDOUT)

        # copy orbfile, if jobdir is not default
        if (self.jobdir != "./"):
            self.orbfile = f'{self.jobdir}{xyzname}.orb'
            os.replace(f'{xyzname}.orb', self.orbfile)
        else:
            self.orbfile = f'{xyzname}.orb'

        self.mol.getEnergies(self.orbfile)

    # Reads energies from an existing HDF5 orb file

 
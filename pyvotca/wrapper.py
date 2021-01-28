"""Votca Wrapper."""
import os
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Dict, Optional

from .molecule import Molecule
from .options import Options
from .xml_editor import edit_xml

__all__ = ["XTP"]


class XTP:

    def __init__(self, mol: Molecule, threads: int = 1, jobname: str = 'dftgwbse',
                 options: Optional[Dict[str, Any]] = {}, jobdir: str = './'):
        self.mol = mol
        self.threads = threads
        self.jobname = jobname
        self.jobdir = jobdir
        self.orbfile = ''
        self.options = Options(options)

    def update_options(self):
        """Merge user options with the defaults."""
        # parsing defaults
        votcashare = os.environ.get('VOTCASHARE')
        default_options = f'{votcashare}/xtp/xml/dftgwbse.xml'
        options = ET.parse(default_options)
        root = options.getroot()
        edit_xml(root, self.options)

        # write out xml
        options.write('dftgwbse.xml')

    def run(self):
        """Just runs xtp_tools with command line call."""
        # update and write the options
        self.update_options()

        # write the XYZfile
        xyzname = self.mol.name
        xyzfile = xyzname + ".xyz"
        self.mol.writeXYZfile(xyzfile)

        """ Runs VOTCA and moves results a job folder, if requested """
        if not Path(self.jobdir).exists():
            os.makedirs(self.jobdir)

        votcacmd = f"xtp_tools -e dftgwbse -o dftgwbse.xml -n {xyzname} -t {self.threads} > {self.jobdir}{self.jobname}.log"

        subprocess.check_output(votcacmd, shell=True, stderr=subprocess.STDOUT)

        # copy orbfile, if jobdir is not default
        if (self.jobdir != "./"):
            self.orbfile = f'{self.jobdir}{self.jobname}.orb'
            os.replace(f'{xyzname}.orb', self.orbfile)
        else:
            self.orbfile = f'{xyzname}.orb'

        self.mol.readORB(self.orbfile)

    # Reads energies from an existing HDF5 orb file

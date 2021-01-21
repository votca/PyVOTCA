from PyVOTCA import XTP
from PyVOTCA import Molecule
from PyVOTCA import Visualization
import numpy as np


# define a molecule
mol=Molecule()

### make it by hand
mol.add_atom("C", 0,0,0)
mol.add_atom("O", 1.2,0,0)

### or read it from existing file
#mol.readXYZfile('CO.xyz')

### get a XTP object
votca=XTP()
#### this allows to change all options
#votca.options['functional'] = 'PBE'
#votca.options['basisset']='cc-pvtz'

# run for the molecule
votca.run(mol)

## only needed, if no run was performed but an existing HDF5 is read
#votca.getEnergies('./CO.orb')

# plotting QP corrections
viz=Visualization(votca)
viz.plotQPcorrections()

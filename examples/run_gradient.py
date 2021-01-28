rom PyVOTCA import XTP
from PyVOTCA import Molecule
from PyVOTCA import NumericalGradient
import numpy as np


# define a molecule
mol=Molecule()

### make it by hand
mol.add_atom("C", 0,0,0)
mol.add_atom("O", 1.3,0.0,0.0)

### get a XTP object
votca=XTP(mol)
#### this allows to change all options
#votca.options['functional'] = 'PBE'
votca.options['basisset']='def2-svp'
votca.options['auxbasisset']='aux-def2-svp'

# run for the molecule it its geometry
votca.run()

# calculate a DFT-GWBSE gradient at the geometry
grad=NumericalGradient(votca,dr=0.001)
grad.getGradient('BSE_singlet',0)

print(mol.getGradient())

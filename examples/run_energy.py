from pyvotca import XTP
from pyvotca import Molecule
from pyvotca import Visualization


# define a molecule
mol = Molecule()

# make it by hand
mol.add_atom("C", 0, 0, 0)
mol.add_atom("O", 1.2, 0, 0)

# or read it from existing file
# mol.readXYZfile('CO.xyz')

# get a XTP object
votca = XTP(mol)
# change basis sets to a smaller one
votca.options['basisset']='def2-svp'
votca.options['auxbasisset']='aux-def2-svp'

# run for the molecule
# votca.run(mol)

# only needed, if no run was performed but an existing HDF5 is read
votca.getEnergies('example.orb')

# Getting the plotting functions
viz = Visualization(votca)
# plotting QP corrections
# viz.plotQPcorrections()
# plotting absorption spectrum
viz.plotAbsorptionGaussian()

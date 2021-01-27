from PyVOTCA import XTP
from PyVOTCA import Molecule
from PyVOTCA import Visualization


# define a molecule
mol = Molecule()

# make it by hand
mol.add_atom("C", 0, 0, 0)
mol.add_atom("O", 1.2, 0, 0)

# or read it from existing file
# mol.readXYZfile('CO.xyz')

# get a XTP object
votca = XTP(mol)
# this allows to change all options
#votca.options['functional'] = 'PBE'
# votca.options['basisset']='cc-pvtz'
votca.options.dftpackage.package.name = "orca"
votca.options.dftpackage.package.executable = "Path/to/Orca"
votca.options.gwbse_engine.gwbse_options.gwbse.mode = 'G0W0'


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

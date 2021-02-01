from pyvotca import Molecule
from pyvotca import Orca
from pyvotca import Electronphonon


# define a molecule
mol=Molecule()

# load xyz 
mol.read_xyz_file('./NPB.xyz')
orca=Orca(mol)

# read gradient from orca
orca.read_gradient('./NPB+.engrad')

# read Hessian from orca
orca.read_hessian('./NPB.hess')

# calculate el-ph couplings and plot
ep = Electronphonon()
ep.calculate_electron_phonon_couplings(mol.elements,mol.hessian,mol.gradient,plot=True)




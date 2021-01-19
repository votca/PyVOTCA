from .numerical_gradient import NumericalGradient as ng
import numpy as np
# INPUT
dr = 0.001  # delta position in Angstrom
# Run multiple simulations in VOTCA with different field directions
ng.run_permut(dr, name='CO', threads=4)
# Get the gradient from all the simulations in Hrtr/Bohr
grad = ng.GradientCalculator(dr,2, "CO")
print(grad.getGradient('dft_tot', energyLevel=7)) # kind of particle/excitation (choices: BSE_singlet, BSE_triplet, QPdiag, QPpert and dft_tot)

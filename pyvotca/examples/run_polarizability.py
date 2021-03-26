#!/usr/bin/env python
"""Example to perform a polarizability calculation."""
import numpy as np

from pyvotca import DFTGWBSE, Molecule, NumericalPolarizability


def run_polarizability():
  mol = Molecule()

  # define a molecule
  mol.add_atom("C", 0, 0, 0)
  mol.add_atom("O", 1.3, 0.0, 0.0)

  # initialize a dftgwbse object
  votca = DFTGWBSE(mol, threads=6)

  # setup the polarizability calculator
  pol = NumericalPolarizability(votca, dE=0.002)
  # run all simulations with different directions of the field
  #pol.run_permut()
  # perform the actual polarizability calculation
  pol_tensor = pol.calc_polarizability("dft_tot")

  print(pol_tensor)


if __name__ == "__main__":
  run_polarizability()
  

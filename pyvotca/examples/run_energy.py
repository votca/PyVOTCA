#!/usr/bin/env python
from pyvotca import DFTGWBSE
from pyvotca import Molecule
from pyvotca import Visualization


def run_energy(save_figure: bool = False):
    """Run energy workflow."""
    # define a molecule
    mol = Molecule()

    # make it by hand
    mol.add_atom("C", 0, 0, 0)
    mol.add_atom("O", 1.2, 0, 0)

    # or read it from existing file
    # mol.readXYZfile('CO.xyz')

    # get a DFTGWBSE object
    dft = DFTGWBSE(mol)
    # change basis sets to a smaller one
    dft.options['basisset'] = 'def2-svp'
    dft.options['auxbasisset'] = 'aux-def2-svp'

    # run for the molecule
    # dft.run(mol)

    # only needed, if no run was performed but an existing HDF5 is read
    dft.mol.read_orb('pyvotca/examples/example.orb')

    # Getting the plotting functions
    viz = Visualization(dft.mol, save_figure=save_figure)
    # plotting QP corrections
    # viz.plot_qp_corrections()
    # plotting absorption spectrum
    viz.plot_absorption_gaussian()


if __name__ == "__main":
    run_energy()

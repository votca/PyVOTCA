"""Entry point to the xtp gradient."""

import argparse
from pathlib import Path

import sys

from .numerical_gradient import NumericalGradient
from .schemas import validate_input


def exists(input_file: str) -> Path:
    """Check if the input file exists."""
    path = Path(input_file)
    if not path.exists():
        raise argparse.ArgumentTypeError(f"{input_file} doesn't exist!")

    return path


def xtp_gradient(args: argparse.Namespace):
    """Compute gradient."""
    # INPUT
    name = args.name
    threads = args.threads
    dr = 0.001  # delta position in Angstrom
    # Run multiple simulations in VOTCA with different field directions
    ng = NumericalGradient(dr, 2, name=name)
    ng.run_permut(dr, 'CO', threads=threads)
    # Get the gradient from all the simulations in Hrtr/Bohr
    grad = ng.GradientCalculator(dr, 2, "CO")
    # kind of particle/excitation (choices: BSE_singlet, BSE_triplet, QPdiag, QPpert and dft_tot)
    print(grad.getGradient('dft_tot', energyLevel=7))


def parse_user_arguments() -> argparse.Namespace:
    """Read the user arguments."""
    parser = argparse.ArgumentParser("xtp_gradient")
    parser.add_argument(
        "-i", "--input", help="Input file in YAML format", type=exists)

    # Read the arguments
    args = parser.parse_args()

    if args.input is None:
        parser.print_help()
        sys.exit()

    return args.input


def main():
    inp = parse_user_arguments()
    args = validate_input(inp)
    print(args)
    # xtp_gradient(args)


if __name__ == "__main__":
    main()

"""Check Orca parsing functionality."""
import numpy as np

from pyvotca.parsers.orca_parsers import (parse_frequencies, parse_gradient,
                                          parse_hessian,
                                          parse_molecular_orbitals,
                                          parse_normal_modes)

from .utils import PATH_TEST


def test_parse_gradient():
    """Check that the gradient is read properly."""
    gradient_file = PATH_TEST / "orca_output" / "CO.engrad"
    grad = parse_gradient(gradient_file)

    # Check that there are not undefined values
    assert not np.all(np.isnan(grad))


def test_parse_hessian():
    """Check that the hessian is read correctly."""
    hessian_file = PATH_TEST / "orca_output" / "methanol.hess"
    hess = parse_hessian(hessian_file)

    # Check that there are not undefined values
    assert not np.all(np.isnan(hess))

    # Check that the hessian is symmetric
    assert np.allclose(hess, hess.T)


def test_read_frequencies():
    """Check that frequencies are read correctly."""
    hessian_file = PATH_TEST / "orca_output" / "methanol.hess"

    freqs = parse_frequencies(hessian_file)

    # Check that there are not undefined values
    assert not np.all(np.isnan(freqs))


def test_read_normal_modes():
    """Check that frequencies are read correctly."""
    hessian_file = PATH_TEST / "orca_output" / "methanol.hess"

    modes = parse_normal_modes(hessian_file)

    # Check that there are not undefined values
    assert not np.all(np.isnan(modes))


def test_read_molecular():
    """Check that molecular orbitals are read correctly."""
    log_file = PATH_TEST / "orca_output" / "methanol.out"

    energies, coefficients = parse_molecular_orbitals(log_file)

    # Check that there are not undefined values
    assert not np.all(np.isnan(energies))
    assert not np.all(np.isnan(coefficients))

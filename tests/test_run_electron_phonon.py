"""Check that the phonon electron example runs properly."""
import numpy as np

from pyvotca.examples.electron_phonon.run_electron_phonon import run_electron_phonon
from .utils import PATH_TEST


def test_run_phonon_electron():
    """Check the phonon electron example."""
    freq, couplings = run_electron_phonon(show=False)

    # Check that the are not NaN values
    assert not np.all(np.isnan(freq))
    assert not np.all(np.isnan(couplings))

    # Check that the results are equal to their references
    ref_freq = np.loadtxt(PATH_TEST / "phonon_electron_frequencies.txt")
    ref_coups = np.loadtxt(PATH_TEST / "phonon_electron_couplings.txt")
    assert np.allclose(ref_freq, freq)
    assert np.allclose(ref_coups, couplings)

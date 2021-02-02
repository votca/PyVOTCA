"""Check that the phonon electron example runs properly."""
import numpy as np

from pyvotca.examples.electron_phonon.run_electron_phonon import run_electron_phonon


def test_run_phonon_electron():
    """Check the phonon electron example."""
    freq, couplings = run_electron_phonon(show=False)

    # Check that the are not NaN values
    assert not np.all(np.isnan(freq))
    assert not np.all(np.isnan(couplings))

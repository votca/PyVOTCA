"""Check that the run_energies examples runs properly."""
from pyvotca.examples.run_energy import run_energy
from .utils import remove_files


def test_run_energy():
    """Check energy run."""
    try:
        run_energy(save_figure=True)
    finally:
        remove_files("*png")

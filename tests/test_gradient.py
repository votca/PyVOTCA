from votcapytools.numerical_gradient import NumericalGradient


def test_instantiation():
    """Instantiate the gradient."""
    NumericalGradient(0.001, 2, "H20")
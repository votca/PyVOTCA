"""Utilities and constants."""
from scipy.constants import physical_constants


class Utils:
    def __init__(self):
        self.h2ev = physical_constants['Hartree energy in eV'][0]

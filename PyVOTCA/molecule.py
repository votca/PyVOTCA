"""Molecule representation."""
import numpy as np
from typing import Union
from pathlib import Path

Pathlike = Union[Path, str]


class Molecule:
    """Molecule definition."""

    def __init__(self):
        self.elements = None
        self.coordinates = None
        self.name = "molecule"

    def add_atom(self, element: str, x: float, y: float, z: float):
        self.elements.append(element)
        self.coordinates.append(np.array([x, y, z]))

    def print(self):
        for (element, coordinates) in zip(self.elements, self.coordinates):
            print(element, coordinates)

    def readXYZfile(self, filename: Pathlike):
        with open(filename, 'r') as handler:
            lines = handler.readlines()

        self.name = Path(filename).stem
        arr = [(row[0], np.array(row[1:], dtype=float)) for row in [
            x.split() for x in lines[2:]]]
        self.elements, self.coordinates = tuple(zip(*arr))

    def writeXYZfile(self, filename: Pathlike):
        """Write the molecule in XYZ format."""
        atoms = "\n".join(f"{elem} {xyz[0]:.4f} {xyz[1]:.4f} {xyz[2]:.4f}" for elem, xyz in zip(
            self.elements, self.coordinates))
        mol = f"""{len(self.elements)}
{self.name} created by PyVOTCA writer
{atoms}
"""
        with open(filename, "w") as xyzfile:
            xyzfile.write(mol)

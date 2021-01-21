
"""Module with the schemas to validate the user input."""

from multiprocessing import cpu_count
from numbers import Integral
from typing import Any, Dict

import yaml
from schema import Optional, Or, Schema, SchemaError

__all__ = ["validate_input"]


input_schema = Schema({
    # Path to the molecule in xyz format
    "molecule": str,

    # Number of Threads to run the application
    Optional("threads", default=cpu_count()): Integral,

    # Functional
    Optional("functional", default="PBE"): str,

    # Basisset
    Optional("basis", default=None): Or(str, None),

    # AuxBasisset
    Optional("auxbasis", default=None): Or(str, None),

    # GW
    Optional("gw", default=None): Or(str, None),

    # BSE
    Optional("bse", default=None): Or(str, None),
})


def validate_input(file_input: str) -> Dict[str, Any]:
    """Check the input validation against an schema."""
    with open(file_input, 'r') as handler:
        dict_input = yaml.load(handler.read(), Loader=yaml.FullLoader)
    try:
        inp = input_schema.validate(dict_input)
        return inp
    except SchemaError as err:
        msg = f"There was an error in the input yaml provided:\n{err}"
        print(msg)
        raise

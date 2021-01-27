"""Module defining the Option class."""
from typing import Any, Dict


class Options(dict):
    """Extend the base class dictionary with a '.' notation.

    example:
    .. code-block:: python
       d = Options({'a': 1})
       d['a'] # 1
       d.a    # 1
       d.b = 3
       d["b"] == 3  # True
    """

    def __init__(self, *args, **kwargs):
        """ Create a recursive Options object"""
        super().__init__(*args, **kwargs)
        for k, v in self.items():
            if isinstance(v, dict):
                self[k] = Options(v)

    def __getattr__(self, attr: str):
        """ Allow `obj.key` notation"""
        return self.get(attr)

    def __setattr__(self, key: str, value: Any):
        """ Allow `obj.key = new_value` notation"""
        self.__setitem__(key, value)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to a normal dictionary."""
        def converter(var):
            return var.to_dict() if isinstance(var, Options) else var

        return {k: converter(v) for k, v in self.items()}
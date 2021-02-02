
from pathlib import Path

import pkg_resources as pkg

__all__ = ["PATH_PYVOTCA", "PATH_TEST", "remove_files"]

# Environment data
PATH_PYVOTCA = Path(pkg.resource_filename('pyvotca', ''))
ROOT = PATH_PYVOTCA.parent

PATH_TEST = ROOT / "tests" / "files"


def remove_files(pattern: str, folder: str = ".") -> None:
    """Remove files with given pattern."""
    path = Path(folder)
    for file in path.glob(pattern):
        file.unlink()

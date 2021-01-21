from pathlib import Path

import pkg_resources as pkg

__all__ = ["PATH_PYVOTCA", "PATH_TEST"]

# Environment data
PATH_PYVOTCA = Path(pkg.resource_filename('PyVOTCA', ''))
ROOT = PATH_PYVOTCA.parent

PATH_TEST = ROOT / "tests" / "files"

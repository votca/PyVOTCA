import copy
from pyvotca.options import Options

data = {'a': 3, 'c': {'d': 42}}


def test_opts():
    """Test the Options class."""
    opts = Options(data)

    assert all((opts.a == data['a'], opts.c['d'] == data['c']['d']))

    # insertion
    opts.elem = 42

    assert opts['elem'] == 42

    # Check copy
    other = copy.deepcopy(opts)
    other.c.d = 23
    assert other.c.d != opts.c.d


def test_undefined_keyword():
    """Check that keywords are created on the fly."""
    opts = Options({"root": 1})
    opts.data.something.foo.bar = 1000

    assert opts.data.something.foo.bar == 1000

    assert opts.data["something"]["foo"]["bar"] == 1000

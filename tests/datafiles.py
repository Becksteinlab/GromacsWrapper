# utilities for data file managements for the tests
"""\
gromacs.tests.datafiles
=======================

In the test code, access a data file "fixtures.dat" in the ``data`` directory with::

  from gromacs.tests.datafiles import datafile

  test_something():
     filepath = datafile("fixtures.dat")
     contents = open(filepath).read()

Basically, wheneever you need the path to the file, wrap the filename in ``datafile()``.

"""

import atexit
from contextlib import ExitStack
import importlib_resources


def datafile(name):
    return importlib_resources.files("tests.data") / name
    file_manager = ExitStack()
    atexit.register(file_manager.close)
    ref = importlib_resources.files("tests.data") / name
    path = file_manager.enter_context(importlib_resources.as_file(ref))
    return path

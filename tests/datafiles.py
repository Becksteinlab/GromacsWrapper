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


import sys
import os.path
from importlib import resources


def datafile(name):
    return (
        resources.files(sys.modules[__name__].__package__)
        .joinpath("data")
        .joinpath(name)
    )

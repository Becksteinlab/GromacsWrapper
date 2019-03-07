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


import os.path
from pkg_resources import resource_filename

def datafile(name):
    return resource_filename(__name__, os.path.join("data", name))



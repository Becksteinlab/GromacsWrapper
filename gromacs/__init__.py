# $Id$
"""\
Gromacs shell
=============

A thin shell around the Gromacs tools for light-weight integration into python
scripts or interactive use in ``ipython``.


Notes
-----

In an ideal world we would have library bindings provided (eg generated via
SWIG wrappers) so that we could directly load the Gromacs analysis code as a
dynamic library. The analysis tools are actually properly compartmentalized
already but the wrappers have not been written. Thus one still has to run
commands using the shell, though now in the disguise of the system() call (or
in python using the os.subprocess
(http://docs.python.org/library/subprocess.html) module in preference to
os.system() (http://docs.python.org/library/os.html)).

To make the integration seamless one should have classes named like the gromacs
programs that take arguments like the programs themselves. In this way the user
will only have to interface with the classes and one can easily switch the
underlying implementation. In addition one can imagine making the analysis tool
MD-package agnostic and switching the 'backend' depending on what kind of
trajectories/topologies one is analyzing.

"""

__all__ = ['tools']

from tools import *

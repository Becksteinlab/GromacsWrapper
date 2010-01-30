# edPDB -- python snippets to edit pdb files
# Copyright (c) 2010 Oliver Beckstein <orbeckst@gmail.com>
#
# Published under the same licence as BioPython.

"""
:mod:`edPDB` -- editing PDB files
=================================

A collection of python snippets to quickly edit pdb files. This is
typically used for setting up system for MD simulations. Of course,
one could use any number of more elegant or powerful tools that
this...

Modules
-------

:mod:`edPDB.cbook`
    Cook-book with short functions that show how to implement basic
    functionality.

:mod:`edPDB.xpdb`
    Extensions to the Bio.PDB class.
"""

try:
    import Bio.PDB as PDB
except ImportError:
    raise ImportError("Biopython's Bio.PDB is absolutely essential. Please install it.")

import logging
# NOTE: logging is still iffy; when I reload I add a new logger each
# time and output is repeated for each reload. Probably should heed
# the advice on logging and libraries in
# http://docs.python.org/library/logging.html?#configuring-logging-for-a-library
class NullHandler(logging.Handler):
    def emit(self, record):
        pass

# this clashes conceptually with the console stuff above: really the
# above needs to be done in application code; in that case the following
# would be enabled:
#
h = NullHandler()
logging.getLogger("edPDB").addHandler(h)
del h

# add standard logging
import log
logger = log.create()


import xpdb
import cbook


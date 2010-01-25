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

import xpdb
import cbook


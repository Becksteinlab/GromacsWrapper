# GromacsWrapper: test_example.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

from __future__ import division, absolute_import, print_function

import os.path

import gromacs.utilities


def test_which(name="cat"):
    path = gromacs.utilities.which(name)
    assert(os.path.basename(gromacs.utilities.which("cat")), name)



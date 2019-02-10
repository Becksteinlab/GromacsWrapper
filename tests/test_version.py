# GromacsWrapper: test_example.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

import pytest

import gromacs

def test_version():
    release = gromacs.__version__
    assert isinstance(release, str)

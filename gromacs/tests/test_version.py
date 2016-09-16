# GromacsWrapper: test_example.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

from __future__ import division, absolute_import, print_function

import pytest

import gromacs.version

def test_get_version():
    release = gromacs.version.get_version()
    assert isinstance(release, str)

def test_get_version_tuple(version=gromacs.version.VERSION):
    release = gromacs.version.get_version_tuple()
    n = len(version)
    assert isinstance(release, tuple)
    assert release[:n] == tuple(map(str, version))

def test_devsuffix(devsuffix="-dev"):
    release = gromacs.version.get_version_tuple()
    # XOR http://stackoverflow.com/questions/432842/how-do-you-get-the-logical-xor-of-two-variables-in-python
    assert gromacs.version.RELEASE is not release[-1].endswith(devsuffix), \
        "{0} should not be present if RELEASE == False: RELEASE = {1}, release = {2}".format(
            devsuffix, gromacs.version.RELEASE, release)

# GromacsWrapper: test_example.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

import gromacs
import pytest
from numpy.testing import assert_equal

from ..datafiles import datafile


@pytest.fixture(
    params=['original', 'nofilename', 'written']
)
def SIMPLE_NDX(request, tmpdir):
    ndx = gromacs.fileformats.ndx.NDX(datafile('simple.ndx'))

    if request.param == 'written':
        out = str(tmpdir.join('out.ndx'))

        ndx.write(out)
        ndx = gromacs.fileformats.ndx.NDX(out)
    elif request.param == 'nofilename':
        ndx = gromacs.fileformats.ndx.NDX()
        ndx.read(datafile('simple.ndx'))

    return ndx


def test_read(SIMPLE_NDX):
    ndx = SIMPLE_NDX

    assert_equal(ndx['Oxygen'], [1, 4, 7])
    assert_equal(ndx['Hydrogen'], [2, 3, 5, 6, 8, 9])


def test_get(SIMPLE_NDX):
    assert_equal(SIMPLE_NDX.get('Oxygen'), [1, 4, 7])


def test_set(SIMPLE_NDX):
    SIMPLE_NDX['Nitrogen'] = [10, 11, 12]

    assert_equal(SIMPLE_NDX['Nitrogen'], [10, 11, 12])


def test_size(SIMPLE_NDX):
    assert len(SIMPLE_NDX) == 2


def test_sizes(SIMPLE_NDX):
    assert SIMPLE_NDX.sizes == {'Oxygen': 3, 'Hydrogen': 6}


def test_groups(SIMPLE_NDX):
    assert list(SIMPLE_NDX.groups) == ['Oxygen', 'Hydrogen']

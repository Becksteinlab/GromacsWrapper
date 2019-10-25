from __future__ import division, absolute_import, print_function

import pytest
from numpy.testing import assert_almost_equal, assert_equal

import numpy as np

import gromacs.fileformats.xpm
from gromacs.formats import XPM

from ..datafiles import datafile

@pytest.fixture
def ssfile():
    return datafile("fileformats/ss.xpm.bz2")

@pytest.fixture
def xpm(ssfile):
    return XPM(filename=ssfile)

@pytest.fixture
def xpm_df(xpm):
    return xpm.to_df()

class TestXPM(object):
    def _run_tests(self, x):
        assert_equal(x.array.shape, (500, 769))
        assert_equal(x.xvalues, np.arange(0, 500))
        assert_equal(x.yvalues, np.arange(1, 770))
        assert_equal((x.array == 'A-Helix').sum(), 292829)

    def test_constructor(self, xpm):
        self._run_tests(xpm)

    def test_read(self, ssfile):
        x = XPM()
        x.read(ssfile)
        self._run_tests(x)

    def test_reversed_by_default(self, xpm):
        assert xpm.reverse

    def test_to_pd(self, xpm_df):
        assert_equal(xpm_df.shape, (500, 770))

    def test_to_pd_types(self, xpm_df):
        time = xpm_df['Time']
        assert len(time) == 500
        assert time.dtype == np.dtype("int64")  # true for this file

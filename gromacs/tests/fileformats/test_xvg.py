from __future__ import division, absolute_import, print_function

import numpy as np

import pytest
from numpy.testing import (assert_almost_equal, assert_array_equal,
                           assert_array_almost_equal, assert_equal, )

from gromacs.formats import XVG


class TestXVG_array():
    def setup(self):
        data = np.random.normal(loc=2.1, scale=0.5, size=(10000, 6))
        data[:, 0] = np.arange(len(data))
        self.x = XVG(array=data.copy(), names="t,a,b,c,d,e")
        self.data = data

    def test_names(self):
        assert_array_equal(self.x.names,  ['t', 'a', 'b', 'c', 'd', 'e'])

    def test_array(self):
        assert_array_equal(self.x.array.shape, self.data.shape)
        assert_array_almost_equal(self.x.array, self.data)

    def test_props(self):
        assert_array_almost_equal(self.x.mean, self.data[1:].mean(axis=1))
        assert_array_almost_equal(self.x.max, self.data[1:].max(axis=1))
        assert_array_almost_equal(self.x.min, self.data[1:].min(axis=1))





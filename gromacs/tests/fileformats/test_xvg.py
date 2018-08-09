from __future__ import division, absolute_import, print_function

import numpy as np
import matplotlib

import pytest
from numpy.testing import assert_almost_equal, assert_equal

import gromacs.fileformats.xvg
from gromacs.formats import XVG

@pytest.fixture(scope="module")
def data():
    data = np.random.normal(loc=2.1, scale=0.5, size=(6, 1000))
    data[0] = 0.1 * np.arange(data.shape[1])
    return data

@pytest.fixture(scope="module")
def correldata(omega=2*np.pi/5):
    t = np.linspace(-100, 100, 10000)
    Y1 = 1 * np.sin(omega*t) - 4
    Y2 = -0.5 * np.sin(0.3*omega*t)
    DY = np.random.normal(scale=2, size=(2, len(t)))
    return np.vstack([t, Y1 + DY[0], Y2 + DY[1]])

class TestXVG_array(object):
    @pytest.fixture
    def xvg(self, data):
        return XVG(array=data.copy(), names="t,a,b,c,d,e")

    def test_names(self, xvg):
        assert_equal(xvg.names,  ['t', 'a', 'b', 'c', 'd', 'e'])

    def test_array(self, xvg, data):
        assert_equal(xvg.array.shape, data.shape)
        assert_almost_equal(xvg.array, data)

    @pytest.mark.parametrize("name",
                             ("mean", "max", "min", "std"))
    def test_props(self, xvg, data, name):
        assert_almost_equal(getattr(xvg, name),
                            getattr(data[1:], name)(axis=1))

    def test_write_read(self, xvg, tmpdir):
        fname = "random.xvg"
        with tmpdir.as_cwd():
            xvg.write(fname)
            newxvg = XVG(filename=fname)
        assert_almost_equal(newxvg.array, xvg.array)
        ## will fail: column names are not written
        # assert_equal(newxvg.names, xvg.names)

    def test_correl(self, correldata):
        xvg = XVG(array=correldata, names="t,y1,y2")
        # FIXME
        # force=True : TypeError: tcorrel() got an unexpected keyword argument 'force'
        xvg.set_correlparameters(nstep=None, ncorrel=25000)

        sigma = xvg.error
        tc = xvg.tc
        assert_equal(sigma.shape, (2,))
        assert_equal(tc.shape, (2,))

    @pytest.mark.parametrize('method', ('mean', 'circmean', 'min', 'max',
                                        'rms', 'percentile', 'smooth',
                                        'error'))
    def test_decimate(self, correldata, method, maxpoints=100):
        xvg = XVG(array=correldata)
        data = xvg.array
        reduced = xvg.decimate(method, data, maxpoints=maxpoints)
        assert_equal(reduced.shape, (len(data), maxpoints))

    def test_plot(self, xvg):
        ax = xvg.plot()
        assert isinstance(ax, matplotlib.axes.Axes)

    def test_errorplot(self, xvg):
        ax = xvg.errorbar(maxpoints=1000)
        assert isinstance(ax, matplotlib.axes.Axes)

    def test_plot_coarsend(self, xvg):
        ax = xvg.plot_coarsened(maxpoints=100)
        assert isinstance(ax, matplotlib.axes.Axes)


def test_break_array():
    angles = np.pi * np.array([-1.9, -1, -1, -0.5, 0, 0.9, 1.5, 2, -2, -1.4])
    expected = np.pi * np.array([-1.9, -1, -1, -0.5, 0, 0.9, 1.5, 2, np.NAN, -2, -1.4])
    other = np.ones_like(angles)
    ma, mother = gromacs.fileformats.xvg.break_array(angles,
                                                     threshold=np.pi, other=other)
    assert isinstance(ma, np.ma.core.MaskedArray)
    assert isinstance(mother, np.ma.core.MaskedArray)
    assert_almost_equal(ma, expected)
    assert len(ma) == len(mother)


# numkit --- time series manipulation and analysis
# Copyright (c) 2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the "Modified BSD Licence" (see COPYING).
"""
:mod:`numkit.timeseries` --- Time series manipulation and analysis
==================================================================

A time series contains of a sequence of time points (typically spaced
equally) and a value for each time point.

.. autoexception:: LowAccuracyWarning


Correlations
------------

.. autofunction:: tcorrel
.. autofunction:: autocorrelation_fft

Autocorrelation time (time when ACF becomes 0 for the first time)::

   R = gromacs.formats.XVG("./md.xvg")
   acf = autocorrelation_fft(R.array[1])
   numpy.where(acf <= 0)[0][0]

Alternatively, fit an exponential to the ACF and extract the time
constant (see :func:`tcorrel`).


Coarse graining time series
---------------------------

The functions in this section are all based on
:func:`regularized_function`. They reduce the number of datapoints in
a time series to *maxpoints* by histogramming the data into
*maxpoints* bins and then applying a function to reduce the data in
each bin. A number of commonly used functions are predefined but it is
straightforward to either use :func:`apply_histogrammed_function` or
:func:`regularized_function` directly. For instance,
:func:`mean_histogrammed_function` is implemented as ::

  def mean_histogrammed_function(t, y, maxpoints):
      return apply_histogrammed_function(numpy.mean, t, y, maxpoints)

More complicated functions can be defined; for instance, one could use
:func:`tcorrel` to compute the correlation time of the data in short
blocks::

   def tc_histogrammed_function(t, y, maxpoints):
      dt = numpy.mean(numpy.diff(t))
      def get_tcorrel(y):
         t = numpy.cumsum(dt*numpy.ones_like(y)) - dt
         results = tcorrel(t, y, nstep=1)
         return results['tc']
      return apply_histogrammed_function(get_tcorrel, t, y, bins=maxpoints)

(This particular function (implemented as :func:`tc_histogrammed_function`) is
not very robust, for instance it has problems when there are only very few data
points in each bin because in this case the auto correlation function is not
well defined.)

.. autofunction:: mean_histogrammed_function
.. autofunction:: rms_histogrammed_function
.. autofunction:: min_histogrammed_function
.. autofunction:: max_histogrammed_function
.. autofunction:: median_histogrammed_function
.. autofunction:: percentile_histogrammed_function
.. autofunction:: error_histogrammed_function
.. autofunction:: circmean_histogrammed_function
.. autofunction:: circstd_histogrammed_function
.. autofunction:: tc_histogrammed_function
.. autofunction:: apply_histogrammed_function
.. autofunction:: regularized_function


Smoothing time series
---------------------

Function :func:`smooth` applies a window kernel to a time series and
smoothes fluctuations. The number of points in the time series stays
the same.

.. autofunction:: smooth
.. autofunction:: smoothing_window_length


"""

from itertools import izip

import numpy
import scipy.signal
import scipy.integrate
import scipy.stats

import warnings
import logging
logger = logging.getLogger("numkit.timeseries")

from numkit import LowAccuracyWarning

def autocorrelation_fft(series, remove_mean=True, paddingcorrection=True,
                        normalize=False, **kwargs):
    """Calculate the auto correlation function.

       autocorrelation_fft(series,remove_mean=False,**kwargs) --> acf

    The time series is correlated with itself across its whole length. Only the
    [0,len(series)[ interval is returned.

    By default, the mean of the series is subtracted and the correlation of the
    fluctuations around the mean are investigated.

    For the default setting remove_mean=True, acf[0] equals the variance of
    the series, acf[0] = Var(series) = <(series - <series>)**2>.

    Optional:

    * The series can be normalized to its 0-th element so that acf[0] == 1.

    * For calculating the acf, 0-padding is used. The ACF should be corrected
      for the 0-padding (the values for larger lags are increased) unless
      mode='valid' is set (see below).

    Note that the series for mode='same'|'full' is inaccurate for long times
    and should probably be truncated at 1/2*len(series)

    :Arguments:
      *series*
        (time) series, a 1D numpy array of length N
      *remove_mean*
        ``False``: use series as is;
        ``True``: subtract mean(series) from series [``True``]
      *paddingcorrection*
        ``False``: corrected for 0-padding; ``True``: return as is it is.
        (the latter is appropriate for periodic signals).
        The correction for element 0=<i<N amounts to a factor N/(N-i). Only
        applied for modes != "valid"       [``True``]
      *normalize*
        ``True`` divides by acf[0] so that the first element is 1;
        ``False`` leaves un-normalized [``False``]
      *mode*
        "full" | "same" | "valid": see :func:`scipy.signal.fftconvolve`
        ["full"]
      *kwargs*
        other keyword arguments for :func:`scipy.signal.fftconvolve`
    """
    kwargs.setdefault('mode','full')

    if len(series.shape) > 2:
        # var/mean below would need proper axis arguments to deal with high dim
        raise TypeError("series must be a 1D array at the moment")

    if remove_mean:
        series = numpy.squeeze(series.astype(float)).copy()   # must copy because de-meaning modifies it
        mean = series.mean()
        series -= mean
    else:
        series = numpy.squeeze(series.astype(float))          # can deal with a view

    ac = scipy.signal.fftconvolve(series,series[::-1,...],**kwargs)

    origin = ac.shape[0]/2        # should work for both odd and even len(series)
    ac = ac[origin:]              # only use second half of the symmetric acf
    assert len(ac) <= len(series), "Oops: len(ac)=%d  len(series)=%d" % (len(ac),len(series))
    if paddingcorrection and  not kwargs['mode'] == 'valid':     # 'valid' was not 0-padded
        # correct for 0 padding
        # XXX: reference? Where did I get this from? (But it makes sense.)
        ac *= len(series)/(len(series) - 1.0*numpy.arange(len(ac)))

    norm = ac[0] or 1.0  # to guard against ACFs of zero arrays
    if not normalize:
        # We use the convention that the ACF is divided by the total time,
        # which makes acf[0] == <series**2> = Var(series) + <series>**2. We do
        # not need to know the time (x) in order to scale the output from the
        # ACF-series accordingly:
        try:
            if remove_mean:
                norm /= numpy.var(series)
            else:
                norm /= numpy.mean(series*series)
        except ZeroDivisionError:
            norm = 1.0
    return ac/norm

def tcorrel(x, y, nstep=100, debug=False):
    """Calculate the correlation time and an estimate of the error of the mean <y>.

    The autocorrelation function f(t) is calculated via FFT on every *nstep* of
    the **fluctuations** of the data around the mean (y-<y>). The normalized
    ACF f(t)/f(0) is assumed to decay exponentially, f(t)/f(0) = exp(-t/tc) and
    the decay constant tc is estimated as the integral of the ACF from the
    start up to its first root.

    See [FrenkelSmit2002]_ `p526`_ for details.

    .. Note:: *nstep* should be set sufficiently large so that there are less
              than ~50,000 entries in the input.

    .. [FrenkelSmit2002] D. Frenkel and B. Smit, Understanding
                         Molecular Simulation. Academic Press, San
                         Diego 2002

    .. _p526: http://books.google.co.uk/books?id=XmyO2oRUg0cC&pg=PA526


    :Arguments:
       *x*
          1D array of abscissa values (typically time)
       *y*
          1D array of the ibservable y(x)
       *nstep*
          only analyze every *nstep* datapoint to speed up calculation
          [100]

    :Returns: dictionary with entries *tc* (decay constant in units of *x*),
              *t0* (value of the first root along x (y(t0) = 0)), *sigma* (error estimate
              for the mean of y, <y>, corrected for correlations in the data).
    """
    if x.shape != y.shape:
        raise TypeError("x and y must be y(x), i.e. same shape")
    _x = x[::nstep]  # do not run acf on all data: takes too long
    _y = y[::nstep]  # and does not improve accuracy
    if len(_y) < 500:  # 500 is a bit arbitrary
        wmsg = "tcorrel(): Only %d datapoints for the chosen nstep=%d; " \
            "ACF will possibly not be accurate." % (len(_y), nstep)
        warnings.warn(wmsg, category=LowAccuracyWarning)
        logger.warn(wmsg)

    acf = autocorrelation_fft(_y, normalize=False)
    try:
        i0 = numpy.where(acf <= 0)[0][0]  # first root of acf
    except IndexError:
        i0 = -1   # use last value as best estimate
    t0 = _x[i0]
    # integral of the _normalized_ acf
    norm = acf[0] or 1.0  # guard against a zero ACF
    tc = scipy.integrate.simps(acf[:i0]/norm, x=_x[:i0])
    # error estimate for the mean [Frenkel & Smit, p526]
    sigma = numpy.sqrt(2*tc*acf[0]/(x[-1] - x[0]))

    result = {'tc':tc, 't0':t0, 'sigma':sigma}
    if debug:
        result['t'] = _x[:i0]
        result['acf'] = acf[:i0]
    return result


def smooth(x, window_len=11, window='flat'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    :Arguments:
        *x*
            the input signal, 1D array

        *window_len*
            the dimension of the smoothing window, always converted to
            an integer (using :func:`int`) and must be odd

        *window*
            the type of window from 'flat', 'hanning', 'hamming',
            'bartlett', 'blackman'; flat window will produce a moving
            average smoothing. If *window* is a :class:`numpy.ndarray` then
            this array is directly used as the window (but it still must
            contain an odd number of points) ["flat"]

    :Returns: the smoothed signal as a 1D array

    :Example:

    Apply a simple moving average to a noisy harmonic signal::

       >>> import numpy as np
       >>> t = np.linspace(-2, 2, 201)
       >>> x = np.sin(t) + np.random.randn(len(t))*0.1
       >>> y = smooth(x)

    .. See Also::

       :func:`numpy.hanning`, :func:`numpy.hamming`,
       :func:`numpy.bartlett`, :func:`numpy.blackman`,
       :func:`numpy.convolve`, :func:`scipy.signal.lfilter`

    Source: based on http://www.scipy.org/Cookbook/SignalSmooth
    """
    windows = {'flat': lambda n: numpy.ones(n, dtype=float),
               'hanning': numpy.hanning,
               'hamming': numpy.hamming,
               'bartlett': numpy.bartlett,
               'blackman': numpy.blackman,
               }
    window_len = int(window_len)

    if isinstance(window, numpy.ndarray):
        window_len = len(window)
        w = numpy.asarray(window, dtype=float)
    else:
        try:
            w = windows[window](window_len)
        except KeyError:
            raise ValueError("Window %r not supported; must be one of %r" %
                             (window, windows.keys()))

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")
    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")
    if window_len % 2 == 0:
        raise ValueError("window_len should be an odd integer")
    if window_len < 3:
        return x

    s = numpy.r_[x[window_len-1:0:-1], x, x[-1:-window_len:-1]]
    y = numpy.convolve(w/w.sum(), s, mode='valid')
    return y[(window_len-1)/2:-(window_len-1)/2]  # take off repeats on ends

def smoothing_window_length(resolution, t):
    """Compute the length of a smooting window of *resolution* time units.

    :Arguments:

       *resolution*
            length in units of the time in which *t* us supplied
       *t*
            array of time points; if not equidistantly spaced, the
            mean spacing is used to compute the window length

    :Returns: odd integer, the size of a window of approximately
              *resolution*

    .. SeeAlso:: :func:`smooth`
    """
    dt = numpy.mean(numpy.diff(t))
    N = int(resolution/dt)
    if N % 2 == 0:
        N += 1
    return N

def mean_histogrammed_function(t, y, **kwargs):
    """Compute mean of data *y* in bins along *t*.

    Returns the mean-regularised function *F* and the centers of the bins.

    .. SeeAlso:: :func:`regularized_function` with *func* = :func:`numpy.mean`
    """
    return apply_histogrammed_function(numpy.mean, t, y, **kwargs)

def rms_histogrammed_function(t, y, **kwargs):
    """Compute root mean square of data *y* in bins along *t*.

    Returns the RMS-regularised function *F* and the centers of the
    bins. *demean* = ``True`` removes the mean first.

    :func:`regularized_function` with *func* = ``sqrt(mean(y*y))``
    """
    def rms(a, demean=kwargs.pop('demean', False)):
        if len(a) == 0:
            return numpy.NAN
        if demean:
            a -= numpy.mean(a)
        return numpy.sqrt(numpy.mean(a*a))
    return apply_histogrammed_function(rms, t, y, **kwargs)

def min_histogrammed_function(t, y, **kwargs):
    """Compute minimum of data *y* in bins along *t*.

    Returns the min-regularised function *F* and the centers of the bins.

    :func:`regularized_function` with *func* = :func:`numpy.min`
    """
    def _min(a):
        if len(a) == 0:
            return numpy.NAN
        return numpy.min(a)
    return apply_histogrammed_function(_min, t, y, **kwargs)

def max_histogrammed_function(t, y, **kwargs):
    """Compute maximum of data *y* in bins along *t*.

    Returns the max-regularised function *F* and the centers of the bins.

    :func:`regularized_function` with *func* = :func:`numpy.max`
    """
    def _max(a):
        if len(a) == 0:
            return numpy.NAN
        return numpy.max(a)
    return apply_histogrammed_function(_max, t, y, **kwargs)

def median_histogrammed_function(t, y, **kwargs):
    """Compute median of data *y* in bins along *t*.

    Returns the median-regularised function *F* and the centers of the bins.

    :func:`regularized_function` with *func* = :func:`numpy.median`
    """
    return apply_histogrammed_function(numpy.median, t, y, **kwargs)

def percentile_histogrammed_function(t, y, **kwargs):
    """Compute the percentile *per* of data *y* in bins along *t*.

    Returns the percentile-regularised function *F* and the centers of
    the bins.

    :Keywords:

      *per*
          percentile as a percentage, e.g. 75 is the value that splits
          the data into the lower 75% and upper 25%; 50 is the median
          [50.0]

      *demean*
          ``True``: remove the mean of the bin data first [``False``]

    :func:`regularized_function` with :func:`scipy.stats.scoreatpercentile`
    """
    def percentile(a, per=kwargs.pop('per', 50.), limit=kwargs.pop('limit', ()),
                   demean=kwargs.pop('demean', False), interpolation_method='fraction'):
        if len(a) == 0:
            return numpy.NAN
        if demean:
            a -= numpy.mean(a)
        return scipy.stats.scoreatpercentile(a, per, limit=limit)
    return apply_histogrammed_function(percentile, t, y, **kwargs)

def tc_histogrammed_function(t, y, **kwargs):
    """Calculate the correlation time in each bin using :func:`tcorrel`.

    .. Warning:: Not well tested and fragile.
    """
    dt = numpy.mean(numpy.diff(t))
    def get_tcorrel(a):
        if len(a) == 0:
            return numpy.NAN
        t = numpy.cumsum(dt*numpy.ones_like(a)) - dt
        results = tcorrel(t, a, nstep=1)
        return results['tc']
    return apply_histogrammed_function(get_tcorrel, t, y, **kwargs)

def error_histogrammed_function(t, y, **kwargs):
    """Calculate the error in each bin using :func:`tcorrel`.

    .. Warning:: Not well tested and fragile.
    """
    dt = numpy.mean(numpy.diff(t))
    def get_tcorrel(a):
        if len(a) == 0:
            return numpy.NAN
        t = numpy.cumsum(dt*numpy.ones_like(a)) - dt
        results = tcorrel(t, a, nstep=1)
        return results['sigma']
    return apply_histogrammed_function(get_tcorrel, t, y, **kwargs)

def circmean_histogrammed_function(t, y, **kwargs):
    """Compute circular mean of data *y* in bins along *t*.

    Returns the circmean-regularised function *F* and the centers of
    the bins.

    *kwargs* are passed to :func:`scipy.stats.morestats.circmean`, in
    particular set the lower bound with *low* and the upper one with
    *high*. The default is [-pi, +pi].

    :func:`regularized_function` with *func* = :func:`scipy.stats.morestats.circmean`

    .. Note:: Data are interpreted as angles in radians.
    """
    low = kwargs.pop('low', -numpy.pi)
    high = kwargs.pop('high', numpy.pi)
    def _circmean(a, low=low, high=high):
        if len(a) == 0:
            return numpy.NAN
        return scipy.stats.morestats.circmean(a, low=low, high=high)
    return apply_histogrammed_function(_circmean, t, y, **kwargs)

def circstd_histogrammed_function(t, y, **kwargs):
    """Compute circular standard deviation of data *y* in bins along *t*.

    Returns the circstd-regularised function *F* and the centers of
    the bins.

    *kwargs* are passed to :func:`scipy.stats.morestats.circmean`, in
    particular set the lower bound with *low* and the upper one with
    *high*. The default is [-pi, +pi].

    :func:`regularized_function` with *func* = :func:`scipy.stats.morestats.circstd`

    .. Note:: Data are interpreted as angles in radians.
    """
    low = kwargs.pop('low', -numpy.pi)
    high = kwargs.pop('high', numpy.pi)
    def _circstd(a, low=low, high=high):
        if len(a) == 0:
            return numpy.NAN
        return scipy.stats.morestats.circstd(a, low=low, high=high)
    return apply_histogrammed_function(_circstd, t, y, **kwargs)

def apply_histogrammed_function(func, t, y, **kwargs):
    """Compute *func* of data *y* in bins along *t*.

    Returns the *func* -regularised function *F(t')* and the centers
    of the bins *t'*.

    .. function:: func(y) -> float

       *func* takes exactly one argument, a numpy 1D array *y* (the
       values in a single bin of the histogram), and reduces it to one
       scalar float.

    """
    F, e = regularized_function(t, y, func, **kwargs)
    return F, 0.5*(e[:-1] + e[1:])

def regularized_function(x, y, func, bins=100, range=None):
    """Compute *func()* over data aggregated in bins.

    ``(x,y) --> (x', func(Y'))``  with ``Y' = {y: y(x) where x in x' bin}``

    First the data is collected in bins x' along x and then *func* is
    applied to all data points Y' that have been collected in the bin.

    .. function:: func(y) -> float

       *func* takes exactly one argument, a numpy 1D array *y* (the
       values in a single bin of the histogram), and reduces it to one
       scalar float.

    .. Note:: *x* and *y* must be 1D arrays.

    :Arguments:
       x
          abscissa values (for binning)
       y
          ordinate values (func is applied)
       func
          a numpy ufunc that takes one argument, func(Y')
       bins
          number or array
       range
          limits (used with number of bins)

    :Returns:
       F,edges
          function and edges (``midpoints = 0.5*(edges[:-1]+edges[1:])``)

    (This function originated as
    :func:`recsql.sqlfunctions.regularized_function`.)
    """
    _x = numpy.asarray(x)
    _y = numpy.asarray(y)

    if len(_x.shape) != 1 or len(_y.shape) != 1:
        raise TypeError("Can only deal with 1D arrays.")

    # setup of bins (taken from numpy.histogram)
    if (range is not None):
        mn, mx = range
        if (mn > mx):
            raise AttributeError('max must be larger than min in range parameter.')

    if not numpy.iterable(bins):
        if range is None:
            range = (_x.min(), _x.max())
        mn, mx = [float(mi) for mi in range]
        if mn == mx:
            mn -= 0.5
            mx += 0.5
        bins = numpy.linspace(mn, mx, bins+1, endpoint=True)
    else:
        bins = numpy.asarray(bins)
        if (numpy.diff(bins) < 0).any():
            raise ValueError('bins must increase monotonically.')

    sorting_index = numpy.argsort(_x)
    sx = _x[sorting_index]
    sy = _y[sorting_index]

    # boundaries in SORTED data that demarcate bins; position in bin_index is the bin number
    bin_index = numpy.r_[sx.searchsorted(bins[:-1], 'left'),
                         sx.searchsorted(bins[-1], 'right')]

    # naive implementation: apply operator to each chunk = sy[start:stop] separately
    #
    # It's not clear to me how one could effectively block this procedure (cf
    # block = 65536 in numpy.histogram) because there does not seem to be a
    # general way to combine the chunks for different blocks, just think of
    # func=median
    F = numpy.zeros(len(bins)-1)  # final function
    F[:] = [func(sy[start:stop]) for start,stop in izip(bin_index[:-1],bin_index[1:])]
    return F,bins

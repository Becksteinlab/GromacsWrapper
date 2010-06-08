# GromacsWrapper -- numerical helper functions for analysis
# Copyright (c) 2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License version 3 (or higher)
"""
:mod:`gromacs.analysis.numkit` --- Mathematical analysis helper functions
=========================================================================

.. autofunction:: Pearson
.. autofunction:: linfit
.. autofunction:: autocorrelation_fft
.. autofunction:: averaged_autocorrelation
.. autofunction:: tcorrel

Examples
--------

Autocorrelation time (time when ACF becomes 0 for the first time)::
  R = gromacs.formats.XVG("./md.xvg")
  acf = mdpow.numkit.autocorrelation_fft(R.array[1])
  where(acf <= 0)[0][0]
 
Alternatively, fit an exponential to the ACF and extract the time constant.

"""

import numpy
import scipy.integrate

import logging
logger = logging.getLogger("gromacs.analysis.numkit")

from gromacs import LowAccuracyWarning

# functions copied from hop.utilities

def Pearson_r(x,y):
    """Pearson's r (correlation coefficient)

      Pearson(x,y) --> correlation coefficient

    *x* and *y* are arrays of same length.
    
    Historical note:
    Naive implementation of Pearson's r:

      Ex = scipy.stats.mean(x)
      Ey = scipy.stats.mean(y)
      covxy = numpy.sum((x-Ex)*(y-Ey))
      r = covxy/math.sqrt(numpy.sum((x-Ex)**2)*numpy.sum((y-Ey)**2))
    """
    return numpy.corrcoef(x,y)[1,0]

def linfit(x,y,dy=[]):
    """Fit a straight line y = a + bx to the data in x and y; errors
    on y should be provided in dy in order to assess the goodness of
    the fit and derive errors on the parameters.

      linfit(x,y[,dy]) --> result_dict

    Fit y = a + bx to the data in x and y by analytically minimizing
    chi^2. dy holds the standard deviations of the individual y_i. If
    dy is not given, they are assumed to be constant (note that in
    this case Q is set to 1 and it is meaningless and chi2 is
    normalised to unit standard deviation on all points!).

    Returns the parameters a and b, their uncertainties sigma_a and
    sigma_b, and their correlation coefficient r_ab; it also returns
    the chi-squared statistic and the goodness-of-fit probability Q
    (that the fit would have chi^2 this large or larger; Q < 10^-2
    indicates that the model is bad --- Q is the probability that a
    value of chi-square as _poor_ as the calculated statistic chi2
    should occur by chance.)
  
    result_dict::
       intercept, sigma_intercept    a +/- sigma_a
       slope, sigma_slope            b +/- sigma_b
       parameter_correlation         correlation coefficient r_ab
                                     between a and b
       chi_square                    chi^2 test statistic
       Q_fit                         goodness-of-fit probability

    Based on 'Numerical Recipes in C', Ch 15.2.
    """
    import scipy.stats

    n = len(x)
    m = len(y)
    if n != m:
        raise ValueError("lengths of x and y must match: %s != %s" % (n, m))
    
    try:
        have_dy = (len(dy) > 0)
    except TypeError:
        have_dy = False

    if not have_dy:
        dy = numpy.ones((n),numpy.float)

    x  = numpy.asarray(x)
    y  = numpy.asarray(y)
    dy = numpy.asarray(dy)

    s2  = dy*dy
    S   = numpy.add.reduce(1/s2)
    Sx  = numpy.add.reduce(x/s2)
    Sy  = numpy.add.reduce(y/s2)
    Sxx = numpy.add.reduce(x*x/s2)
    Sxy = numpy.add.reduce(x*y/s2)

    t   = (x - Sx/S)/dy
    Stt = numpy.add.reduce(t*t)

    b = numpy.add.reduce(t*y/dy)/Stt
    a = (Sy - Sx*b)/S

    sa = numpy.sqrt((1 + (Sx*Sx)/(S*Stt))/S)
    sb = numpy.sqrt(1/Stt)

    covab = -Sx/(S*Stt)
    r = covab/(sa*sb)

    chi2 = numpy.add.reduce(((y-a-b*x)/dy)**2)
    if not have_dy:
        # estimate error if none were provided
        sigmadata = numpy.sqrt(chi2/(n-2))
        sa *= sigmadata
        sb *= sigmadata
        Q = 1.0
    else:
        Q = scipy.stats.chisqprob(chi2,n-2)

    return {"intercept":a,"slope":b,
            "sigma_intercept":sa,"sigma_slope":sb,
            "parameter_correlation":r, "chi_square":chi2, "Q":Q}

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

    * The series is can be normalized to its 0-th element so that acf[0] == 1.

    * For calculating the acf, 0-padding is used. The ACF shoudl be corrected
      for the 0-padding (the values for larger lags are increased) unless
      mode='valid' (see below).

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
    import scipy.signal
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

def tcorrel(x,y,nstep=100,debug=False):
    """Calculate the correlation time and an estimate of the error of <y>.

    The autocorrelation function f(t) is calculated via FFT on every *nstep* of
    the **fluctuations** of the data around the mean (y-<y>). The normalized
    ACF f(t)/f(0) is assumed to decay exponentially, f(t)/f(0) = exp(-t/tc) and
    the decay constant tc is estimated as the integral of the ACF from the
    start up to its first root.

    See Frenkel and Smit, Academic Press, San Diego 2002, p526.

    .. Note:: *nstep* should be set sufficiently large so that there are less
              than ~50,000 entries in the input.

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
        import warnings
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




class FitFunc(object):
    """Fit a function f to data (x,y) using the method of least squares.

    The function is fitted when the object is created, using
    :func:`scipy.optimize.leastsq`. One must derive from the base class
    :class:`FitFunc` and override the :meth:`FitFunc.f_factory` (including
    the definition of an appropriate local :func:`fitfunc` function) and
    :meth:`FitFunc.initial_values` appropriately. See the examples for a
    linear fit :class:`FitLin`, a 1-parameter exponential fit :class:`FitExp`,
    or a 3-parameter double exponential fit :class:`FitExp2`.

    The object provides two attributes
     :attr:`FitFunc.parameters`
           list of parameters of the fit
     :attr:`FitFunc.message`
           message from :func:`scipy.optimize.leastsq`

    After a successful fit, the fitted function can be applied to any data (a
    1D-numpy array) with :meth:`FitFunc.fit`.
     
    """
    def __init__(self,x,y):
        import scipy.optimize
        _x = numpy.asarray(x)
        _y = numpy.asarray(y)
        p0 = self.initial_values()
        fitfunc = self.f_factory()
        def errfunc(p,x,y):
            return  fitfunc(p,x) - y     # residuals        
        p,msg = scipy.optimize.leastsq(errfunc,p0[:],args=(_x,_y))
        try:
            p[0]
            self.parameters = p
        except (TypeError,IndexError,):
            # TypeError for int p, IndexError for numpy scalar (new scipy)
            self.parameters = [p]
        self.message = msg

    def f_factory(self):
        """Stub for fit function factory, which returns the fit function.
        Override for derived classes.
        """
        def fitfunc(p,x):
            # return f(p,x); should be a numpy ufunc
            raise NotImplementedError("base class must be extended for each fit function")
        return fitfunc

    def initial_values(self):
        """List of initital guesses for all parameters p[]"""
        # return [1.0, 2.0, 0.5]
        raise NotImplementedError("base class must be extended for each fit function")    

    def fit(self,x):
        """Applies the fit to all *x* values"""
        fitfunc = self.f_factory()
        return fitfunc(self.parameters,numpy.asarray(x))

class FitExp(FitFunc):
    """y = f(x) = exp(-p[0]*x)"""
    def f_factory(self):
        def fitfunc(p,x):
            return numpy.exp(-p[0]*x)   # exp(-B*x)
        return fitfunc
    def initial_values(self):
        return [1.0]
    def __repr__(self):
        return "<FitExp "+str(self.parameters)+">"

class FitExp2(FitFunc):
    """y = f(x) = p[0]*exp(-p[1]*x) + (1-p[0])*exp(-p[2]*x)"""
    def f_factory(self):
        def fitfunc(p,x):
            return p[0]*numpy.exp(-p[1]*x) + (1-p[0])*numpy.exp(-p[2]*x)
        return fitfunc
    def initial_values(self):
        return [0.5,0.1,1e-4]
    def __repr__(self):
        return "<FitExp2"+str(self.parameters)+">"

class FitLin(FitFunc):
    """y = f(x) = p[0]*x + p[1]"""
    def f_factory(self):
        def fitfunc(p,x):
            return p[0]*x + p[1]
        return fitfunc
    def initial_values(self):
        return [1.0,0.0]
    def __repr__(self):
        return "<FitLin"+str(self.parameters)+">"

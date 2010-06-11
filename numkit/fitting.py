# numkit --- data fitting 
# Copyright (c) 2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the "Modified BSD Licence" (see COPYING).
"""
:mod:`numkit.fitting` --- Fitting data
======================================

.. autofunction:: Pearson_r
.. autofunction:: linfit

.. autoclass:: FitFunc
   :members:
.. autoclass:: FitLin
.. autoclass:: FitExp
.. autoclass:: FitExp2

"""

import numpy
import logging
logger = logging.getLogger("numkit.fitting")


def Pearson_r(x,y):
    """Pearson's r (correlation coefficient).

       Pearson(x,y) --> correlation coefficient

    *x* and *y* are arrays of same length.
    
    Historical note -- Naive implementation of Pearson's r ::
      Ex = scipy.stats.mean(x)
      Ey = scipy.stats.mean(y)
      covxy = numpy.sum((x-Ex)*(y-Ey))
      r = covxy/math.sqrt(numpy.sum((x-Ex)**2)*numpy.sum((y-Ey)**2))
    """
    return numpy.corrcoef(x,y)[1,0]

def linfit(x,y,dy=[]):
    """Fit a straight line y = a + bx to the data in *x* and *y*.

    Errors on y should be provided in dy in order to assess the
    goodness of the fit and derive errors on the parameters.

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
  
    :Returns: result_dict with components

       intercept, sigma_intercept    
           a +/- sigma_a
       slope, sigma_slope
           b +/- sigma_b
       parameter_correlation
           correlation coefficient r_ab between a and b
       chi_square                    
           chi^2 test statistic
       Q_fit
           goodness-of-fit probability

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


from scipy.integrate.quadrature import tupleset

def simps_error(y, x=None, dx=1, axis=-1, even='avg'):
    """Error on integral evaluated with `Simpson's rule`_ from errors of points, y.

    Evaluate the integral with :func:`scipy.integrate.simps`. For a
    given vector *y* of errors on the function values, the error on
    the integral is calculated via propagation of errors.

    .. Note: If the spacing is not equal then the error propagation is
             made with the *approximation* of a constant spacing equal
             to the average *x* spacing. This should be fixed in a future
             release.

    :Arguments:
      *y*
         errors for the tabulated values of the integrand f
      *x*
         values of abscissa at which f was tabulated (can be ``None``
         and then *dx* should be provided)
      *dx*
         constant spacing of the abscissa
      *axis*
         axis in *y* along which the data lies
      *even*
         see :func:`scipy.integrate.simps` ('avg', 'first', 'last')

    .. _Simpson's rule: http://mathworld.wolfram.com/SimpsonsRule.html
    """
    # copied basic structure from scipy.integrate.quadrature.simps
    y = numpy.asarray(y)
    nd = len(y.shape)
    N = y.shape[axis]
    if not x is None:
        x = numpy.asarray(x)
        if len(x.shape) == 1:
            shapex = numpy.ones(nd)
            shapex[axis] = x.shape[0]
            saveshape = x.shape
            returnshape = 1
            x=x.reshape(tuple(shapex))
        elif len(x.shape) != len(y.shape):
            raise ValueError, "If given, shape of x must be 1-d or the " \
                  "same as y."
        if x.shape[axis] != N:
            raise ValueError, "If given, length of x along axis must be the " \
                  "same as y."        
    if N % 2 == 0:
        val = 0.0      # holds trapezoidal error**2
        result = 0.0   # holds Simposon error**2
        slice1 = (slice(None),)*nd
        slice2 = (slice(None),)*nd
        if not even in ['avg', 'last', 'first']:
            raise ValueError, \
                  "Parameter 'even' must be 'avg', 'last', or 'first'."
        # Compute using Simpson's rule on first intervals
        if even in ['avg', 'first']:
            slice1 = tupleset(slice1, axis, -1)
            slice2 = tupleset(slice2, axis, -2)
            if not x is None:
                last_dx = x[slice1] - x[slice2]
            val += (0.5*last_dx)**2 * (y[slice1]+y[slice2])**2  # trapz. error
            result = _simps_error(y,0,N-3,x,dx,axis)
        # Compute using Simpson's rule on last set of intervals
        if even in ['avg', 'last']:
            slice1 = tupleset(slice1, axis, 0)
            slice2 = tupleset(slice2, axis, 1)
            if not x is None:
                first_dx = x[tuple(slice2)] - x[tuple(slice1)]
            val += (0.5*first_dx)**2 * (y[slice2]+y[slice1])**2 # trapz. error
            result += _simps_error(y,1,N-2,x,dx,axis)
        if even == 'avg':
            val /= 2.0**2     # error propagation on y=(a1+a2)/2 gives
            result /= 2.0**2  # dy**2 = (da1**2+da2**2)/2**2 
            # (although not quite correct as errors da1 and da2 are not independent)
        result = result + val
    else:
        result = _simps_error(y,0,N-2,x,dx,axis)
    if returnshape:
        x = x.reshape(saveshape)
    return numpy.sqrt(result)

def _simps_error(y,start,stop,x,dx,axis):
    """Squared error on Simpson's rule integration.

    .. Note: If the spacing is not equal then the error propagation is
             made with the *approximation* of a constant spacing equal
             to the average spacing.

    :Arguments:
       *y*
          errors at function values
       *start*, *stop*
          first and last index at which a Simpson 3-point interval starts
       *x*
          abscissa values (provide if spacing not equal)
       *dx*
          constant spacing (is overridden by *dx*)
       *axis*
          axis in *y* along which the data lie
    """
    nd = len(y.shape)
    if start is None:
        start = 0
    step = 2
    all = (slice(None),)*nd
    slice0 = tupleset(all, axis, slice(start, stop, step))
    slice1 = tupleset(all, axis, slice(start+1, stop+1, step))
    slice2 = tupleset(all, axis, slice(start+2, stop+2, step))

    # Simpson error propgation for points 0 <= i <= 2M
    # error**2 = (h/3)**2 + sum_k=1^M dy[2k]**2 + (4*dy[2k-1])**2 + dy[2k-2]**2

    if x is None:  # Even spaced Simpson's rule.        
        h = dx
    else:
        # simplified assumption: all spacings are equal... or at least so
        # on average
        h = numpy.mean(x[start+1:stop+2] - x[start:stop+1])

    result = numpy.add.reduce((h/3.0)**2 * ((y[slice0])**2+(4*y[slice1])**2+(y[slice2])**2), axis)

    return result

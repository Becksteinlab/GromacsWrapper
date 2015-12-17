# numkit --- data fitting
# Copyright (c) 2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the "Modified BSD Licence" (see COPYING).
"""
:mod:`numkit.fitting` --- Fitting data
======================================

The module contains functions to do least square fits of functions of
one variable f(x) to data points (x,y).

Example
-------

For example, to fit a un-normalized Gaussian with :class:`FitGauss` to
data distributed with mean 5.0 and standard deviation 3.0::

   from numkit.fitting import FitGauss
   import numpy, numpy.random

   # generate suitably noisy data
   mu, sigma = 5.0, 3.0
   Y,edges = numpy.histogram(sigma*numpy.random.randn(10000), bins=100, density=True)
   X = 0.5*(edges[1:]+edges[:-1]) + mu

   g = FitGauss(X, Y)

   print(g.parameters)
   # [ 4.98084541  3.00044102  1.00069061]
   print(numpy.array([mu, sigma, 1]) - g.parameters)
   # [ 0.01915459 -0.00044102 -0.00069061]

   import matplotlib.pyplot as plt
   plt.plot(X, Y, 'ko', label="data")
   plt.plot(X, g.fit(X), 'r-', label="fit")

.. figure:: /numkit/FitGauss.png
   :scale: 40 %
   :alt: Gaussian fit with data points

   A Gaussian (red) was fit to the data points (black circles) with
   the :class:`numkit.fitting.FitGauss` class.

If the initial parameters for the least square optimization do not
lead to a solution then one can provide customized starting values in
the *parameters* keyword argument::

   g = FitGauss(X, Y, parameters=[10, 1, 1])

The *parameters* have different meaning for the different fit
functions; the documentation for each function shows them in the
context of the fit function.


Creating new fit functions
--------------------------

New fit function classes can be derived from :class:`FitFunc`. The
documentation and the methods :meth:`FitFunc.f_factory` and
:meth:`FitFunc.initial_values` must be overriden. For example, the
class :class:`FitGauss` is implemented as ::

   class FitGauss(FitFunc):
       '''y = f(x) = p[2] * 1/sqrt(2*pi*p[1]**2) * exp(-(x-p[0])**2/(2*p[1]**2))'''
       def f_factory(self):
           def fitfunc(p,x):
               return p[2] * 1.0/(p[1]*numpy.sqrt(2*numpy.pi)) * numpy.exp(-(x-p[0])**2/(2*p[1]**2))
           return fitfunc
       def initial_values(self):
           return [0.0,1.0,0.0]

The function to be fitted is defined in :func:`fitfunc`. The
parameters are accessed as ``p[0]``, ``p[1]``, ... For each parameter,
a suitable initial value must be provided.


Functions and classes
---------------------

.. autofunction:: Pearson_r
.. autofunction:: linfit

.. autoclass:: FitFunc
   :members:
.. autoclass:: FitLin
.. autoclass:: FitExp
.. autoclass:: FitExp2
.. autoclass:: FitGauss

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

def linfit(x,y,dy=None):
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
    if dy is None:
        dy = []
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
    def __init__(self,x,y,parameters=None):
        import scipy.optimize
        _x = numpy.asarray(x)
        _y = numpy.asarray(y)
        p0 = self._get_initial_values(parameters)
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

    def _get_initial_values(self, parameters=None):
        p0 = numpy.asarray(self.initial_values())
        if parameters is not None:
            try:
                p0[:] = parameters
            except ValueError:
                raise ValueError("Wrong number of custom initital values %r, should be something like %r" % (parameters, p0))
        return p0

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

class FitGauss(FitFunc):
    """y = f(x) = p[2] * 1/sqrt(2*pi*p[1]**2) * exp(-(x-p[0])**2/(2*p[1]**2))

    Fits an un-normalized gaussian (height scaled with parameter p[2]).

    * p[0] == mean $\mu$
    * p[1] == standard deviation $\sigma$
    * p[2] == scale $a$
    """
    def f_factory(self):
        def fitfunc(p,x):
            return p[2] * 1.0/(p[1]*numpy.sqrt(2*numpy.pi)) * numpy.exp(-(x-p[0])**2/(2*p[1]**2))
        return fitfunc
    def initial_values(self):
        return [0.0,1.0,0.0]
    def __repr__(self):
        return "<FitGauss"+str(self.parameters)+">"

# numkit --- integration
# Copyright (c) 2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the "Modified BSD Licence" (see COPYING).
"""
:mod:`numkit.integration` --- Numerical integration of data
===========================================================

.. SeeAlso:: :mod:`scipy.integrate`

.. autofunction:: simps_error
.. autoexception:: LowAccuracyWarning

"""

import numpy
import scipy.integrate
from scipy.integrate.quadrature import tupleset

import logging
logger = logging.getLogger("numkit.integration")

class LowAccuracyWarning(Warning):
    """Warns that results may possibly have low accuracy."""

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

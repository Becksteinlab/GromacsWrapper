# numkit --- integration
# Copyright (c) 2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the "Modified BSD Licence" (see COPYING).
"""
:mod:`numkit.integration` --- Numerical integration of data
===========================================================

.. SeeAlso:: :mod:`scipy.integrate`

.. autofunction:: simps_error

"""

import numpy
import scipy.integrate
from scipy.integrate.quadrature import tupleset

import logging
logger = logging.getLogger("numkit.integration")

def simps_error(dy, x=None, dx=1, axis=-1, even='avg'):
    """Error on integral evaluated with `Simpson's rule`_ from errors of points, *dy*.

    Evaluate the integral with :func:`scipy.integrate.simps`. For a
    given vector *dy* of errors on the function values, the error on
    the integral is calculated via `propagation of errors`_ from the
    `Newton-Cotes formula`_ for the 3rd `Lagrange interpolating
    polynomial`_. The results are exact for the cases of even spacing
    *dx*; for uneven spacing we currently average all spacings (exact
    solution is in the works...)

    :Arguments:
      *dy*
         errors for the tabulated values of the integrand f
      *x*
         values of abscissa at which f was tabulated (can be ``None``
         and then *dx* should be provided)
      *dx*
         constant spacing of the abscissa
      *axis*
         axis in *dy* along which the data lies
      *even*
         see :func:`scipy.integrate.simps` ('avg', 'first', 'last')

    .. _Simpson's rule: http://mathworld.wolfram.com/SimpsonsRule.html
    .. _propagation of errors: http://mathworld.wolfram.com/ErrorPropagation.html
    .. _`Newton-Cotes formula`: http://mathworld.wolfram.com/Newton-CotesFormulas.html
    .. _Lagrange interpolating polynomial: http://mathworld.wolfram.com/LagrangeInterpolatingPolynomial.html
    """
    # copied basic structure from scipy.integrate.quadrature.simps
    dy = numpy.asarray(dy)
    nd = len(dy.shape)
    N = dy.shape[axis]

    if not x is None:
        x = numpy.asarray(x)
        if len(x.shape) == 1:
            shapex = numpy.ones(nd)
            shapex[axis] = x.shape[0]
            saveshape = x.shape
            returnshape = True
            x=x.reshape(tuple(shapex))
        elif len(x.shape) != len(dy.shape):
            raise ValueError, "If given, shape of x must be 1-d or the " \
                  "same as dy."
        if x.shape[axis] != N:
            raise ValueError, "If given, length of x along axis must be the " \
                  "same as dy."
    else:
        last_dx = first_dx = dx
        avglast_dx = avgfirst_dx = dx
        returnshape = False

    if N % 2 == 0:
        val = 0.0      # holds trapezoidal error**2
        result = 0.0   # holds Simpson error**2
        correction = 0.0  # overlap correction for combining trapz and Simpson
        slice1 = (slice(None),)*nd
        slice2 = (slice(None),)*nd
        if not even in ['avg', 'last', 'first']:
            raise ValueError, \
                  "Parameter 'even' must be 'avg', 'last', or 'first'."
        if even in ['avg', 'first']:
            # Compute using Simpson's rule on first intervals
            slice0 = tupleset(slice1, axis, -3)
            slice1 = tupleset(slice1, axis, -2)
            slice2 = tupleset(slice2, axis, -1)
            if not x is None:
                last_dx =   x[slice2] - x[slice1]
                penult_dx = x[slice1] - x[slice0]
                avglast_dx = 0.5*(penult_dx + last_dx)
            val += _trapz_2pt_error2(dy[slice1], dy[slice2], last_dx)
            result += _simps_error2(dy,0,N-3,x,dx,axis)
            correction += _trapz_simps_overlap2(dy[slice1], avglast_dx)  # avglast_dx not strictly correct for x!=None!
        if even in ['avg', 'last']:
            # Compute using Simpson's rule on last set of intervals
            slice1 = tupleset(slice1, axis, 0)
            slice2 = tupleset(slice2, axis, 1)
            slice3 = tupleset(slice1, axis, 2)
            if not x is None:
                first_dx = x[slice2] - x[slice1]
                second_dx = x[slice3] - x[slice2]
                avgfirst_dx = 0.5*(second_dx + first_dx)
            val += _trapz_2pt_error2(dy[slice1], dy[slice2], first_dx)
            result += _simps_error2(dy,1,N-2,x,dx,axis)
            correction += _trapz_simps_overlap2(dy[slice2], avgfirst_dx)  # avgfirst_dx not strictly correct for x!=None!
        if even == 'avg':
            # simply average the variances of 'odd' and 'even' because the two
            # data sets are not independent; hence it would be wrong to
            # estimate it as  avg**2 = 0.5*sqrt(da1**2+da2**2)
            val /= 2.0
            result /= 2.0
            correction /= 2.0
        result = result + val + correction  # err_simps**2 + err_trapez**2 + correction
    else:
        result = _simps_error2(dy,0,N-2,x,dx,axis)
    if returnshape:
        x = x.reshape(saveshape)
    return numpy.sqrt(result)

def _trapz_simps_overlap2(dy, dx):
    """Correction term in the squared error when combining trapezoidal and Simpson's rule.

    Only exact for equal spacing *dx* left and right of *dy*.

    err^2 = (h/6)^2 ((3 Df0)^2 + ((3+2)Df1)^2 + (8Df2)^2 + (4Df3)^2 + ...)
                      |-- trapz ---| |--------- Simpson ---------------

    (3+2)^2 = 3^2 + 2^2 + 12    <--- 12 is the "overlap" correction
    """
    return (dx/6)**2 * 12 * dy**2

def _simps_error2(dy,start,stop,x,dx,axis):
    """Squared error on Simpson's rule integration.

    :Arguments:
       *dy*
          errors at function values
       *start*, *stop*
          first and last index at which a Simpson 3-point interval starts
       *x*
          abscissa values (provide if spacing not equal)
       *dx*
          constant spacing (is overridden by *dx*)
       *axis*
          axis in *dy* along which the data lie
    """
    nd = len(dy.shape)
    if start is None:
        start = 0
    step = 2
    all = (slice(None),)*nd
    # check that stop is appropriate !!!
    slice2k  = tupleset(all, axis, slice(start, stop, step))     # does NOT include last point stop+1
    slice2k1 = tupleset(all, axis, slice(start+1, stop+1, step)) 
    slice0   = tupleset(all, axis, slice(start, start+1))          # first point
    sliceN   = tupleset(all, axis, slice(stop+1, stop+2))          # last point

    Df2k  = dy[slice2k]    # 2k
    Df2k1 = dy[slice2k1]   # 2k+1
    Df0   = dy[slice0]     # first = 0
    DfN   = dy[sliceN]     # last = N

    if x is None:  # Even spaced Simpson's rule.
        # Simpson error propagation for I = (h/3)*(f_1 + 4f_2 + 2f_3 + 4f_4 + ... + 4f_{N-1} + f_N)
        # error**2 = (h/3)**2 * [sum_k=0^(N-1)/2 ((2*Df_2k)**2 + (4*Df_2k+1)**2) - 3*(Df_0**2 + Df_N**2)]
        #          = (h/3)**2 * [sum_k=0^(N-3)/2 ((2*Df_2k)**2 + (4*Df_2k+1)**2) - 3*Df_0**2 + Df_N**2]
        # with the last term being the correction for the end points
        #
        # In order to work with nd-arrays I need to get the axis of  '- 3*Df0**2 + DfN**2'
        # but my numpy-foo is weak and I don't know how to extract this axis nicely: Hence
        # sum the single term with the axis argument :-(
        result = (dx/3.0)**2 * (numpy.add.reduce(((2*Df2k)**2 + (4*Df2k1)**2), axis) \
                                    + numpy.add.reduce(-3*Df0**2 + DfN**2, axis))
    else:
        logger.warning("Approximating Simpson integration statistical error with the average spacing.")
        dx = numpy.diff(x).mean()
        result = _simps_error2(dy,start,stop,None,dx,axis)

    return result

def _trapz_2pt_error2(dy1, dy2, dx):
    """Squared error on the trapezoidal rule for two points.

    For errors dy1 and dy2 and spacing dx."""
    return (0.5*dx)**2 * (dy1**2 + dy2**2)

def _trapz_error2(dy,start,stop,x,dx):
    # not tested
    nd = len(dy.shape)
    if start is None:
        start = 0
    step = 1
    all = (slice(None),)*nd
    slicek   = tupleset(all, axis, slice(start, stop+2, step))     # all points; stop+2 ??
    slice0   = tupleset(all, axis, slice(start, start+1))          # first point
    sliceN   = tupleset(all, axis, slice(stop+1, stop+2))          # last point

    Dfk  = dy[slicek]    # 2k
    Df0   = dy[slice0]     # first = 0
    DfN   = dy[sliceN]     # last = N

    if x is None:
        # from error propagation of I = (h/2) * (f_1 + 2f_2 + 2f_3 + ... + 2f_{N-1) + f_N)
        result = (0.5*dx)**2 * (4 * numpy.add.reduce(Dfk**2, axis) - 3*(Df0**2 + DfN**2))
    else:
        raise NotImplementedError
    return result



def _naive_simps_error2(dy,start,stop,x,dx,axis):
    """Simple squared error on Simpson's rule integration.

    The implementation assumes that the errors on the individual
    intervals are independent and hence their squares can be
    summed. This is wrong. DO NOT USE.

    :Arguments:
       *dy*
          errors at function values
       *start*, *stop*
          first and last index at which a Simpson 3-point interval starts
       *x*
          abscissa values (provide if spacing not equal)
       *dx*
          constant spacing (is overridden by *dx*)
       *axis*
          axis in *dy* along which the data lie
    """
    import warnings
    warnings.warn("The Simpson error from un-even intervals is not correct. _naive_simps_error2() will be removed.",
                  category=DeprecationWarning)

    nd = len(dy.shape)
    if start is None:
        start = 0
    step = 2
    all = (slice(None),)*nd
    slice0 = tupleset(all, axis, slice(start, stop, step))
    slice1 = tupleset(all, axis, slice(start+1, stop+1, step))
    slice2 = tupleset(all, axis, slice(start+2, stop+2, step))

    Df1 = dy[slice0]   # 2k
    Df2 = dy[slice1]   # 2k-1
    Df3 = dy[slice2]   # 2k-2
    
    if x is None:  # Even spaced Simpson's rule.
        # Simpson error propgation for points 0 <= i <= 2M
        # error**2 = (h/3)**2 * sum_k=1^M dy[2k]**2 + (4*dy[2k-1])**2 + dy[2k-2]**2
        # error**2 = (h/3)**2 * (Df1**2 + 4*Df2**2 + Df**2)
        h = dx
        result = numpy.add.reduce((h/3.0)**2 * (Df1**2+(4*Df2)**2+Df3**2), axis)
    else:
        # This is too naive: adds the squared errors for every individual Simpson-triplet but that misses
        # the correlations in errors from the overlapping endpoints.

        # for spacing h1 and h2 (evaluated in Sage from the integral of the Lagrange interpolating 
        # polynomials and propagation of errors)
        # Z = h1**2 * h2 + h1 * h2**2
        # error**2 = 1/(6*Z)**2 * ((2*h1**3*h2 + 3*h1**2*h2**2 - h2**4)**2*Df_1**2 
        #               + (h1**4 - 3*h1**2*h2**2 - 2*h1*h2**3)**2*Df_3**2 
        #               + (h1**4 + 4*h1**3*h2 + 6*h1**2*h2**2 + 4*h1*h2**3 + h2**4)**2*Df_2**2)
        h1 = x[slice1] - x[slice0]  # check this!
        h2 = x[slice2] - x[slice1]
        #print 'h1', h1.shape, h1
        #print 'h2', h2.shape, h2

        Z = h1**2 * h2 + h1 * h2**2

        #print 'Z ', Z
        #print 'Df1 ', Df1
        #print 'Df2 ', Df2
        #print 'Df2 ', Df2
        
        result = numpy.add.reduce(1/(6*Z)**2 * 
                                  (((2*h1**3*h2 + 3*h1**2*h2**2 - h2**4) * Df1)**2 \
                                       + ((h1**4 + 4*h1**3*h2 + 6*h1**2*h2**2 + 4*h1*h2**3 + h2**4) * Df2)**2 \
                                       + ((h1**4 - 3*h1**2*h2**2 - 2*h1*h2**3) * Df3)**2))
    return result


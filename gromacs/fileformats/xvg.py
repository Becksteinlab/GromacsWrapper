# -*- encoding: utf-8 -*-
# GromacsWrapper: formats.py
# Copyright (c) 2009-2012 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.
"""
Simple xmgrace XVG file format
==============================

Gromacs produces graphs in the `xmgrace`_ ("xvg") format. These are
simple multi-column data files. The class :class:`XVG` encapsulates
access to such files and adds a number of methods to access the data
(as NumPy arrays), compute aggregates, or quickly plot it.

.. _xmgrace: http://plasma-gate.weizmann.ac.il/Grace/

The :class:`XVG` class is useful beyond reading xvg files. With the
*array* keyword or the :meth:`XVG.set` method one can load data from
an array instead of a file. The array should be simple "NXY" data
(typically: first column time or position, further columns scalar
observables). The data should be a NumPy :class:`numpy.ndarray` array
``a`` with :attr:`~numpy.ndarray.shape` ``(M, N)`` where *M*-1 is the
number of observables and *N* the number of observations, e.g.the
number of time points in a time series. ``a[0]`` is the time or
position and ``a[1:]`` the *M*-1 data columns.


Errors
------

The :attr:`XVG.error` attribute contains the statistical error for
each timeseries. It is computed from the standard deviation of the
fluctuations from the mean and their correlation time. The parameters
for the calculations of the correlation time are set with
:meth:`XVG.set_correlparameters`.

.. SeeAlso:: :func:`numkit.timeseries.tcorrel`


Plotting
--------

The :meth:`XVG.plot` and :meth:`XVG.errorbar` methods are set up to
produce graphs of multiple columns simultaneously. It is
typically assumed that the first column in the selected (sub)array
contains the abscissa ("x-axis") of the graph and all further columns
are plotted against the first one.

Data selection
~~~~~~~~~~~~~~

Plotting from :class:`XVG` is fairly flexible as one can always pass
the *columns* keyword to select which columns are to be
plotted. Assuming that the data contains ``[t, X1, X2, X3]``, then one
can

1) plot all observable columns (X1 to X3) against t::

     xvg.plot()

2) plot only X2 against t::

     xvg.plot(columns=[0,2])

3) plot X2 and X3 against t::

     xvg.plot(columns=[0,2,3])

4) plot X1 against X3::

     xvg.plot(columns=[2,3])


Coarse grainining of data
~~~~~~~~~~~~~~~~~~~~~~~~~

It is also possible to *coarse grain the data* for plotting (which
typically results in visually smoothing the graph because noise is
averaged out).

Currently, two alternative algorithms to produce "coarse grained"
(decimated) graphs are implemented and can be selected with the
*method* keyword for the plotting functions in conjuction with
*maxpoints* (the number of points to be plotted):

1) **mean** histogram (default) --- bin the data (using
   :func:`numkit.timeseries.regularized_function` and compute the
   mean for each bin. Gives the exact number of desired points
   but the time data are whatever the middle of the bin is.

2) **smooth** subsampled --- smooth the data with a running average
   (other windows like Hamming are also possible) and then pick data
   points at a stepsize compatible with the number of data points
   required. Will give exact times but not the exact number of data
   points.

For simple test data, both approaches give very similar output.

For the special case of periodic data such as angles, one can use the
circular mean ("circmean") to coarse grain. In this case, jumps across
the -180º/+180º boundary are added as masked datapoints and no line is
drawn across the jump in the plot. (Only works with the simple
:meth:`XVG.plot` method at the moment, errorbars or range plots are
not implemented yet.)

.. SeeAlso:: :meth:`XVG.decimate`

Examples
--------

In this example we generate a noisy time series of a sine wave. We
store the time, the value, and an error. (In a real example, the
value might be the mean over multiple observations and the error might
be the estimated error of the mean.)

  >>> import numpy as np
  >>> import gromacs.formats
  >>> X = np.linspace(-10,10,50000)
  >>> yerr = np.random.randn(len(X))*0.05
  >>> data = np.vstack((X, np.sin(X) + yerr, np.random.randn(len(X))*0.05))
  >>> xvg = gromacs.formats.XVG(array=data)

Plot value for *all* time points::

  >>> xvg.plot(columns=[0,1], maxpoints=None, color="black", alpha=0.2)

Plot bin-averaged (decimated) data with the errors, over 1000 points::

  >>> xvg.errorbar(maxpoints=1000, color="red")

(see output in Figure :ref:`Plot of Raw vs Decimated data <figure-xvg-decimated-label>`)

.. _figure-xvg-decimated-label:

.. figure:: xvg_decimated.*
   :figwidth: 40%
   :scale: 70%
   :alt: plot of a raw noisy sin(x) graph versus its decimated version
   :align: right

   **Plot of Raw vs Decimated data.** Example of plotting raw data
   (sine on 50,000 points, gray) versus the decimated graph (reduced
   to 1000 points, red line). The errors were also decimated and
   reduced to the errors within the 5% and the 95% percentile. The
   decimation is carried out by histogramming the data in the desired
   number of bins and then the data in each bin is reduced by either
   :func:`numpy.mean` (for the value) or
   :func:`scipy.stats.scoreatpercentile` (for errors).

In principle it is possible to use other functions to decimate the
data. For :meth:`XVG.plot`, the *method* keyword can be changed (see
:meth:`XVG.decimate` for allowed *method* values). For
:meth:`XVG.errorbar`, the method to reduce the data values (typically
column 1) is fixed to "mean" but the errors (typically columns 2 and 3)
can also be reduced with *error_method* = "rms".

If one wants to show the variation of the raw data together with the
decimated and smoothed data then one can plot the percentiles of the
deviation from the mean in each bin:

   >>> xvg.errorbar(columns=[0,1,1], maxpoints=1000, color="blue", demean=True)

The *demean* keyword indicates if fluctuations from the mean are
regularised [#demean]_. The method :meth:`XVG.plot_coarsened`
automates this approach and can plot coarsened data selected by the
*columns* keyword.

.. rubric:: Footnotes

.. [#demean] When *error_method* = "percentile" is selected for
             :meth:`XVG.errorbar` then *demean* does not actually
             force a regularisation of the fluctuations from the
             mean. Instead, the (symmetric) percentiles are computed
             on the full data and the error ranges for plotting are
             directly set to the percentiles. In this way one can
             easily plot the e.g. 10th-percentile to 90th-percentile
             band (using keyword *percentile* = 10).



Classes and functions
---------------------

.. autoclass:: XVG
   :members:

.. autofunction:: break_array

"""


from __future__ import with_statement

import os, errno
import re
import warnings
from collections import OrderedDict as odict

import numpy

from gromacs.exceptions import (ParseError, MissingDataError,
                                MissingDataWarning, AutoCorrectionWarning)
import gromacs.utilities as utilities


import numkit.timeseries

import logging


class XVG(utilities.FileUtils):
    """Class that represents the numerical data in a grace xvg file.

    All data must be numerical. :const:`NAN` and :const:`INF` values are
    supported via python's :func:`float` builtin function.

    The :attr:`~XVG.array` attribute can be used to access the the
    array once it has been read and parsed. The :attr:`~XVG.ma`
    attribute is a numpy masked array (good for plotting).

    Conceptually, the file on disk and the XVG instance are considered the same
    data. Whenever the filename for I/O (:meth:`XVG.read` and :meth:`XVG.write`) is
    changed then the filename associated with the instance is also changed to reflect
    the association between file and instance.

    With the *permissive* = ``True`` flag one can instruct the file reader to skip
    unparseable lines. In this case the line numbers of the skipped lines are stored
    in :attr:`XVG.corrupted_lineno`.

    A number of attributes are defined to give quick access to simple statistics such as

     - :attr:`~XVG.mean`: mean of all data columns
     - :attr:`~XVG.std`: standard deviation
     - :attr:`~XVG.min`: minimum of data
     - :attr:`~XVG.max`: maximum of data
     - :attr:`~XVG.error`: error on the mean, taking correlation times into
       account (see also :meth:`XVG.set_correlparameters`)
     - :attr:`~XVG.tc`: correlation time of the data (assuming a simple
       exponential decay of the fluctuations around the mean)

    These attributes are numpy arrays that correspond to the data columns,
    i.e. :attr:`XVG.array`[1:].

    .. Note::

       - Only simple XY or NXY files are currently supported, *not*
         Grace files that contain multiple data sets separated by '&'.
       - Any kind of formatting (i.e. :program:`xmgrace` commands) is discarded.
    """

    #: Default extension of XVG files.
    default_extension = "xvg"

    #: Aim for plotting around that many points
    maxpoints_default = 10000

    # logger: for pickling to work, this *must* be class-level and
    # cannot be done in __init__() (because we cannot pickle self.logger)
    logger = logging.getLogger('gromacs.formats.XVG')

    #: If :attr:`XVG.savedata` is ``False`` then any attributes in
    #: :attr:`XVG.__pickle_excluded` are *not* pickled as they are but simply
    #: pickled with the default value.
    __pickle_excluded = {'__array': None}   # note class name un-mangling in __getstate__()!

    #: Default color cycle for :meth:`XVG.plot_coarsened`:
    #: ``['black', 'red', 'blue', 'orange', 'magenta', 'cyan', 'yellow', 'brown', 'green']``
    default_color_cycle = ['black', 'red', 'blue', 'orange', 'magenta', 'cyan', 'yellow', 'brown', 'green']

    def __init__(self, filename=None, names=None, array=None, permissive=False, **kwargs):
        """Initialize the class from a xvg file.

        :Arguments:
              *filename*
                    is the xvg file; it can only be of type XY or
                    NXY. If it is supplied then it is read and parsed
                    when :attr:`XVG.array` is accessed.
              *names*
                    optional labels for the columns (currently only
                    written as comments to file); string with columns
                    separated by commas or a list of strings
              *array*
                    read data from *array* (see :meth:`XVG.set`)
              *permissive*
                    ``False`` raises a :exc:`ValueError` and logs and errior
                    when encountering data lines that it cannot parse.
                    ``True`` ignores those lines and logs a warning---this is
                    a risk because it might read a corrupted input file [``False``]
              *stride*
                    Only read every *stride* line of data [1].
              *savedata*
                    ``True`` includes the data (:attr:`XVG.array`` and
                    associated caches) when the instance is pickled (see
                    :mod:`pickle`); this is oftens not desirable because the
                    data are already on disk (the xvg file *filename*) and the
                    resulting pickle file can become very big. ``False`` omits
                    those data from a pickle. [``False``]
              *metadata*
                    dictionary of metadata, which is not touched by the class

        """
        self.__array = None           # cache for array (BIG) (used by XVG.array)
        self.__cache = {}             # cache for computed results
        self.savedata = kwargs.pop('savedata', False)
        if filename is not None:
            self._init_filename(filename)  # note: reading data from file is delayed until required
        if names is None:
            self.names = []
        else:
            try:
                self.names = names.split(',')
            except AttributeError:
                self.names = names
        self.metadata = kwargs.pop('metadata', {})  # reserved for user data
        self.permissive = permissive
        self.stride = kwargs.pop('stride', 1)
        self.corrupted_lineno = None      # must parse() first before this makes sense
        # default number of data points for calculating correlation times via FFT
        self.ncorrel = kwargs.pop('ncorrel', 25000)
        self.__correlkwargs = {}          # see set_correlparameters()

        if array is not None:
            self.set(array)

    def read(self, filename=None):
        """Read and parse xvg file *filename*."""
        self._init_filename(filename)
        self.parse()

    def write(self, filename=None):
        """Write array to xvg file *filename* in NXY format.

        .. Note:: Only plain files working at the moment, not compressed.
        """
        self._init_filename(filename)
        with utilities.openany(self.real_filename, 'w') as xvg:
            xvg.write("# xmgrace compatible NXY data file\n"
                      "# Written by gromacs.formats.XVG()\n")
            xvg.write("# :columns: %r\n" % self.names)
            for xyy in self.array.T:
                xyy.tofile(xvg, sep=" ", format="%-8s")  # quick and dirty ascii output...--no compression!
                xvg.write('\n')

    @property
    def array(self):
        """Represent xvg data as a (cached) numpy array.

        The array is returned with column-first indexing, i.e. for a data file with
        columns X Y1 Y2 Y3 ... the array a will be a[0] = X, a[1] = Y1, ... .
        """
        if self.__array is None:
            self.parse()
        return self.__array

    @property
    def ma(self):
        """Represent data as a masked array.

        The array is returned with column-first indexing, i.e. for a data file with
        columns X Y1 Y2 Y3 ... the array a will be a[0] = X, a[1] = Y1, ... .

        inf and nan are filtered via :func:`numpy.isfinite`.
        """
        a = self.array
        return numpy.ma.MaskedArray(a, mask=numpy.logical_not(numpy.isfinite(a)))

    @property
    def mean(self):
        """Mean value of all data columns."""
        return self.array[1:].mean(axis=1)

    @property
    def std(self):
        """Standard deviation from the mean of all data columns."""
        return self.array[1:].std(axis=1)

    @property
    def min(self):
        """Minimum of the data columns."""
        return self.array[1:].min(axis=1)

    @property
    def max(self):
        """Maximum of the data columns."""
        return self.array[1:].max(axis=1)

    def _tcorrel(self, nstep=100, **kwargs):
        """Correlation "time" of data.

        The 0-th column of the data is interpreted as a time and the
        decay of the data is computed from the autocorrelation
        function (using FFT).

        .. SeeAlso:: :func:`numkit.timeseries.tcorrel`
        """
        from gromacs.analysis.collections import Collection
        t = self.array[0,::nstep]
        r = Collection([numkit.timeseries.tcorrel(t, Y, nstep=1, **kwargs) for Y in self.array[1:,::nstep]])
        return r

    def set_correlparameters(self, **kwargs):
        """Set and change the parameters for calculations with  correlation functions.

        The parameters persist until explicitly changed.

        :Keywords:
           *nstep*
               only process every *nstep* data point to speed up the FFT; if
               left empty a default is chosen that produces roughly 25,000 data
               points (or whatever is set in *ncorrel*)
           *ncorrel*
               If no *nstep* is supplied, aim at using *ncorrel* data points for
               the FFT; sets :attr:`XVG.ncorrel` [25000]
           *force*
               force recalculating correlation data even if cached values are
               available
           *kwargs*
               see :func:`numkit.timeseries.tcorrel` for other options

        .. SeeAlso: :attr:`XVG.error` for details and references.
        """
        self.ncorrel = kwargs.pop('ncorrel', self.ncorrel) or 25000
        nstep = kwargs.pop('nstep', None)
        if nstep is None:
            # good step size leads to ~25,000 data points
            nstep = len(self.array[0])/float(self.ncorrel)
            nstep = int(numpy.ceil(nstep))  # catch small data sets
        kwargs['nstep'] = nstep
        self.__correlkwargs.update(kwargs)  # only contains legal kw for numkit.timeseries.tcorrel or force
        return self.__correlkwargs

    def _correlprop(self, key, **kwargs):
        kwargs = self.set_correlparameters(**kwargs)
        if not self.__cache.get('tcorrel',None) or kwargs.pop('force', False):
            self.__cache['tcorrel'] = self._tcorrel(**kwargs)
        return numpy.array(self.__cache['tcorrel'].get(key).tolist())

    @property
    def error(self):
        """Error on the mean of the data, taking the correlation time into account.

        See [FrenkelSmit2002]_ `p526`_:

           error = sqrt(2*tc*acf[0]/T)

        where acf() is the autocorrelation function of the fluctuations around
        the mean, y-<y>, tc is the correlation time, and T the total length of
        the simulation.

        .. [FrenkelSmit2002] D. Frenkel and B. Smit, Understanding
                             Molecular Simulation. Academic Press, San
                             Diego 2002

        .. _p526: http://books.google.co.uk/books?id=XmyO2oRUg0cC&pg=PA526
        """
        return self._correlprop('sigma')

    @property
    def tc(self):
        """Correlation time of the data.

        See :meth:`XVG.error` for details.
        """
        return self._correlprop('tc')

    def parse(self, stride=None):
        """Read and cache the file as a numpy array.

        Store every *stride* line of data; if ``None`` then the class default is used.

        The array is returned with column-first indexing, i.e. for a data file with
        columns X Y1 Y2 Y3 ... the array a will be a[0] = X, a[1] = Y1, ... .
        """
        if stride is None:
            stride = self.stride
        self.corrupted_lineno = []
        irow  = 0  # count rows of data
        # cannot use numpy.loadtxt() because xvg can have two types of 'comment' lines
        with utilities.openany(self.real_filename) as xvg:
            rows = []
            ncol = None
            for lineno,line in enumerate(xvg):
                line = line.strip()
                if len(line) == 0:
                    continue
                if "label" in line and "xaxis" in line:
                        self.xaxis = line.split('"')[-2]
                if "label" in line and "yaxis" in line:
                        self.yaxis = line.split('"')[-2]
                if line.startswith("@ legend"):
                                        if not "legend" in self.metadata: self.metadata["legend"] = []
                                        self.metadata["legend"].append(line.split("legend ")[-1])
                if line.startswith("@ s") and "subtitle" not in line:
                                        name = line.split("legend ")[-1].replace('"','').strip()
                                        self.names.append(name)
                if line.startswith(('#', '@')) :
                                        continue
                if line.startswith('&'):
                    raise NotImplementedError('%s: Multi-data not supported, only simple NXY format.'
                                              % self.real_filename)
                # parse line as floats
                try:
                    row = map(float, line.split())
                except:
                    if self.permissive:
                        self.logger.warn("%s: SKIPPING unparsable line %d: %r",
                                         self.real_filename, lineno+1, line)
                        self.corrupted_lineno.append(lineno+1)
                        continue
                    self.logger.error("%s: Cannot parse line %d: %r",
                                      self.real_filename, lineno+1, line)
                    raise
                # check for same number of columns as in previous step
                if not ncol is None and len(row) != ncol:
                    if self.permissive:
                        self.logger.warn("%s: SKIPPING line %d with wrong number of columns: %r",
                                         self.real_filename, lineno+1, line)
                        self.corrupted_lineno.append(lineno+1)
                        continue
                    errmsg = "%s: Wrong number of columns in line %d: %r" % (self.real_filename, lineno+1, line)
                    self.logger.error(errmsg)
                    raise IOError(errno.ENODATA, errmsg, self.real_filename)
                # finally: a good line
                if irow % stride == 0:
                    ncol = len(row)
                    rows.append(row)
                irow += 1
        try:
            self.__array = numpy.array(rows).transpose()    # cache result
        except:
            self.logger.error("%s: Failed reading XVG file, possibly data corrupted. "
                              "Check the last line of the file...", self.real_filename)
            raise
        finally:
            del rows     # try to clean up as well as possible as it can be massively big

    def to_df(self):
        import pandas as _pd
        return _pd.DataFrame(self.array.T, columns=[self.xaxis] + (self.names if len(self.names) else [self.yaxis]) , dtype=float)

    def set(self, a):
        """Set the *array* data from *a* (i.e. completely replace).

        No sanity checks at the moment...
        """
        self.__array = numpy.asarray(a)

    def plot(self, **kwargs):
        """Plot xvg file data.

        The first column of the data is always taken as the abscissa
        X. Additional columns are plotted as ordinates Y1, Y2, ...

        In the special case that there is only a single column then this column
        is plotted against the index, i.e. (N, Y).

        :Keywords:
          *columns* : list
               Select the columns of the data to be plotted; the list
               is used as a numpy.array extended slice. The default is
               to use all columns. Columns are selected *after* a transform.
          *transform* : function
               function ``transform(array) -> array`` which transforms
               the original array; must return a 2D numpy array of
               shape [X, Y1, Y2, ...] where X, Y1, ... are column
               vectors.  By default the transformation is the
               identity [``lambda x: x``].
          *maxpoints* : int
               limit the total number of data points; matplotlib has issues processing
               png files with >100,000 points and pdfs take forever to display. Set to
               ``None`` if really all data should be displayed. At the moment we simply
               decimate the data at regular intervals. [10000]
          *method*
               method to decimate the data to *maxpoints*, see :meth:`XVG.decimate`
               for details
          *color*
               single color (used for all plots); sequence of colors
               (will be repeated as necessary); or a matplotlib
               colormap (e.g. "jet", see :mod:`matplotlib.cm`). The
               default is to use the :attr:`XVG.default_color_cycle`.
          *kwargs*
               All other keyword arguments are passed on to :func:`pylab.plot`.
        """
        from itertools import izip, cycle
        import matplotlib.cm, matplotlib.colors
        import pylab

        columns = kwargs.pop('columns', Ellipsis)         # slice for everything
        maxpoints = kwargs.pop('maxpoints', self.maxpoints_default)
        transform = kwargs.pop('transform', lambda x: x)  # default is identity transformation
        method = kwargs.pop('method', "mean")

        if columns is Ellipsis or columns is None:
            columns = numpy.arange(self.array.shape[0])
        if len(columns) == 0:
            raise MissingDataError("plot() needs at least one column of data")

        if len(self.array.shape) == 1 or self.array.shape[0] == 1:
            # special case: plot against index; plot would do this automatically but
            # we'll just produce our own xdata and pretend that this was X all along
            a = numpy.ravel(self.array)
            X = numpy.arange(len(a))
            a = numpy.vstack((X, a))
            columns = [0] + [c+1 for c in columns]
        else:
            a = self.array

        color = kwargs.pop('color', self.default_color_cycle)
        try:
            cmap = matplotlib.cm.get_cmap(color)
            colors = cmap(matplotlib.colors.Normalize()(numpy.arange(len(columns[1:]), dtype=float)))
        except TypeError:
            colors = cycle(utilities.asiterable(color))

        # (decimate/smooth o slice o transform)(array)
        a = self.decimate(method, numpy.asarray(transform(a))[columns], maxpoints=maxpoints)

        # now deal with infs, nans etc AFTER all transformations (needed for plotting across inf/nan)
        ma = numpy.ma.MaskedArray(a, mask=numpy.logical_not(numpy.isfinite(a)))

        # finally plot (each column separately to catch empty sets)
        for column, color in izip(xrange(1,len(columns)), colors):
            if len(ma[column]) == 0:
                warnings.warn("No data to plot for column %(column)d" % vars(), category=MissingDataWarning)
            kwargs['color'] = color
            pylab.plot(ma[0], ma[column], **kwargs)   # plot all other columns in parallel

    def plot_coarsened(self, **kwargs):
        """Plot data like :meth:`XVG.plot` with the range of **all** data shown.

        Data are reduced to *maxpoints* (good results are obtained
        with low values such as 100) and the actual range of observed
        data is plotted as a translucent error band around the mean.

        Each column in *columns* (except the abscissa, i.e. the first
        column) is decimated (with :meth:`XVG.decimate`) and the range
        of data is plotted alongside the mean using
        :meth:`XVG.errorbar` (see for arguments). Additional
        arguments:

        :Kewords:
           *maxpoints*
                number of points (bins) to coarsen over
           *color*
                single color (used for all plots); sequence of colors
                (will be repeated as necessary); or a matplotlib
                colormap (e.g. "jet", see :mod:`matplotlib.cm`). The
                default is to use the :attr:`XVG.default_color_cycle`.
           *method*
                Method to coarsen the data. See :meth:`XVG.decimate`

        The *demean* keyword has no effect as it is required to be ``True``.

        .. SeeAlso:: :meth:`XVG.plot`, :meth:`XVG.errorbar` and :meth:`XVG.decimate`
        """
        from itertools import izip, cycle
        import matplotlib.cm, matplotlib.colors

        columns = kwargs.pop('columns', Ellipsis)         # slice for everything
        if columns is Ellipsis or columns is None:
            columns = numpy.arange(self.array.shape[0])
        if len(columns) < 2:
            raise MissingDataError("plot_coarsened() assumes that there is at least one column "
                                   "of data for the abscissa and one or more for the ordinate.")

        color = kwargs.pop('color', self.default_color_cycle)
        try:
            cmap = matplotlib.cm.get_cmap(color)
            colors = cmap(matplotlib.colors.Normalize()(numpy.arange(len(columns[1:]), dtype=float)))
        except TypeError:
            colors = cycle(utilities.asiterable(color))

        t = columns[0]
        kwargs['demean'] = True
        for column, color in izip(columns[1:], colors):
            kwargs['color'] = color
            self.errorbar(columns=[t, column, column], **kwargs)

    def errorbar(self, **kwargs):
        """errorbar plot for a single time series with errors.

        Set *columns* keyword to select [x, y, dy] or [x, y, dx, dy],
        e.g. ``columns=[0,1,2]``. See :meth:`XVG.plot` for
        details. Only a single timeseries can be plotted and the user
        needs to select the appropriate columns with the *columns*
        keyword.

        By default, the data are decimated (see :meth:`XVG.plot`) for
        the default of *maxpoints* = 10000 by averaging data in
        *maxpoints* bins.

        x,y,dx,dy data can plotted with error bars in the x- and
        y-dimension (use *filled* = ``False``).

        For x,y,dy use *filled* = ``True`` to fill the region between
        y±dy. *fill_alpha* determines the transparency of the fill
        color. *filled* = ``False`` will draw lines for the error
        bars. Additional keywords are passed to
        :func:`pylab.errorbar`.

        By default, the errors are decimated by plotting the 5% and
        95% percentile of the data in each bin. The percentile can be
        changed with the *percentile* keyword; e.g. *percentile* = 1
        will plot the 1% and 99% perentile (as will *percentile* =
        99).

        The *error_method* keyword can be used to compute errors as
        the root mean square sum (*error_method* = "rms") across each
        bin instead of percentiles ("percentile"). The value of the
        keyword *demean* is applied to the decimation of error data
        alone.

        .. SeeAlso::

           :meth:`XVG.plot` lists keywords common to both methods.
        """
        import pylab

        color = kwargs.pop('color', 'black')
        filled = kwargs.pop('filled', True)
        fill_alpha = kwargs.pop('fill_alpha', 0.2)

        kwargs.setdefault('capsize', 0)
        kwargs.setdefault('elinewidth', 1)
        kwargs.setdefault('ecolor', color)
        kwargs.setdefault('alpha', 0.3)
        kwargs.setdefault('fmt', None)

        columns = kwargs.pop('columns', Ellipsis)         # slice for everything
        maxpoints = kwargs.pop('maxpoints', self.maxpoints_default)
        transform = kwargs.pop('transform', lambda x: x)  # default is identity transformation
        method = kwargs.pop('method', "mean")
        if method != "mean":
            raise NotImplementedError("For errors only method == 'mean' is supported.")
        error_method = kwargs.pop('error_method', "percentile")  # can also use 'rms' and 'error'
        percentile = numpy.abs(kwargs.pop('percentile', 95.))
        demean = kwargs.pop('demean', False)

        # order: (decimate/smooth o slice o transform)(array)
        try:
            data = numpy.asarray(transform(self.array))[columns]
        except IndexError:
            raise MissingDataError("columns %r are not suitable to index the transformed array, possibly not eneough data" % columns)
        if data.shape[-1] == 0:
            raise MissingDataError("There is no data to be plotted.")
        a = numpy.zeros((data.shape[0], maxpoints), dtype=numpy.float64)
        a[0:2] = self.decimate("mean", data[0:2], maxpoints=maxpoints)
        error_data = numpy.vstack((data[0], data[2:]))
        if error_method == "percentile":
            if percentile > 50:
                upper_per = percentile
                lower_per = 100 - percentile
            else:
                upper_per = 100 - percentile
                lower_per = percentile
            # demean generally does not make sense with the percentiles (but for analysing
            # the regularised data itself we use this as a flag --- see below!)
            upper = a[2:] = self.decimate("percentile", error_data, maxpoints=maxpoints,
                                          per=upper_per, demean=False)[1:]
            lower = self.decimate("percentile", error_data, maxpoints=maxpoints,
                                  per=lower_per, demean=False)[1:]
        else:
            a[2:] = self.decimate(error_method, error_data, maxpoints=maxpoints, demean=demean)[1:]
            lower = None

        # now deal with infs, nans etc AFTER all transformations (needed for plotting across inf/nan)
        ma = numpy.ma.MaskedArray(a, mask=numpy.logical_not(numpy.isfinite(a)))
        if lower is not None:
            mlower = numpy.ma.MaskedArray(lower, mask=numpy.logical_not(numpy.isfinite(lower)))

        # finally plot
        X = ma[0]          # abscissa set separately
        Y = ma[1]
        try:
            kwargs['yerr'] = ma[3]
            kwargs['xerr'] = ma[2]
        except IndexError:
            try:
                kwargs['yerr'] = ma[2]
            except IndexError:
                raise TypeError("Either too few columns selected or data does not have a error column")

        if filled:
            # can only plot dy
            if error_method == "percentile":
                if demean:
                    # signal that we are looking at percentiles of an observable and not error
                    y1 = mlower[-1]
                    y2 = kwargs['yerr']
                else:
                    # percentiles of real errors (>0)
                    y1 = Y - mlower[-1]
                    y2 = Y + kwargs['yerr']
            else:
                y1 = Y - kwargs['yerr']
                y2 = Y + kwargs['yerr']
            pylab.fill_between(X, y1, y2, color=color, alpha=fill_alpha)
        else:
            if error_method == "percentile":
                # errorbars extend to different lengths;
                if demean:
                    kwargs['yerr'] = numpy.vstack((mlower[-1], kwargs['yerr']))
                else:
                    kwargs['yerr'] = numpy.vstack((Y - mlower[-1], Y + kwargs['yerr']))
                try:
                    # xerr only makes sense when the data is a real
                    # error so we don't even bother with demean=?
                    kwargs['xerr'] = numpy.vstack((X - mlower[0], X + kwargs['xerr']))
                except (KeyError, IndexError):
                    pass
            pylab.errorbar(X, Y, **kwargs)

        # clean up args for plot
        for kw in "yerr", "xerr", "capsize", "ecolor", "elinewidth", "fmt":
            kwargs.pop(kw, None)
        kwargs['alpha'] = 1.0

        pylab.plot(X, Y, color=color, **kwargs)

    def decimate(self, method, a, **kwargs):
        """Decimate data *a* to *maxpoints* using *method*.

        If *a* is a 1D array then it is promoted to a (2, N) array
        where the first column simply contains the index.

        If the array contains fewer than *maxpoints* points or if
        *maxpoints* is ``None`` then it is returned as it is. The
        default for *maxpoints* is 10000.

        Valid values for the reduction *method*:

          * "mean", uses :meth:`XVG.decimate_mean` to coarse grain by
            averaging the data in bins along the time axis

          * "circmean", uses :meth:`XVG.decimate_circmean` to coarse
            grain by calculating the circular mean of the data in bins
            along the time axis. Use additional keywords *low* and
            *high* to set the limits. Assumes that the data are in
            degrees.

          * "min" and "max* select the extremum in each bin

          * "rms", uses :meth:`XVG.decimate_rms` to coarse grain by
            computing the root mean square sum of the data in bins
            along the time axis (for averaging standard deviations and
            errors)

          * "percentile" with keyword *per*: :meth:`XVG.decimate_percentile`
            reduces data in each bin to the percentile *per*

          * "smooth", uses :meth:`XVG.decimate_smooth` to subsample
            from a smoothed function (generated with a running average
            of the coarse graining step size derived from the original
            number of data points and *maxpoints*)

        :Returns: numpy array ``(M', N')`` from a ``(M', N)`` array
                  with ``M' == M`` (except when ``M == 1``, see above)
                  and ``N' <= N`` (``N'`` is *maxpoints*).
        """
        methods = {'mean': self.decimate_mean,
                   'min': self.decimate_min,
                   'max': self.decimate_max,
                   'smooth': self.decimate_smooth,
                   'rms': self.decimate_rms,
                   'percentile': self.decimate_percentile,
                   'error': self.decimate_error,  # undocumented, not working well
                   'circmean': self.decimate_circmean,
                   }
        if len(a.shape) == 1:
            # add first column as index
            # (probably should do this in class/init anyway...)
            X = numpy.arange(len(a))
            a = numpy.vstack([X, a])
        ny = a.shape[-1]   # assume 1D or 2D array with last dimension varying fastest
        maxpoints = kwargs.pop("maxpoints", 10000)
        if maxpoints is None or ny <= maxpoints:
            return a
        return methods[method](a, maxpoints, **kwargs)

    def decimate_mean(self, a, maxpoints, **kwargs):
        """Return data *a* mean-decimated on *maxpoints*.

        Histograms each column into *maxpoints* bins and calculates
        the weighted average in each bin as the decimated data, using
        :func:`numkit.timeseries.mean_histogrammed_function`. The coarse grained
        time in the first column contains the centers of the histogram
        time.

        If *a* contains <= *maxpoints* then *a* is simply returned;
        otherwise a new array of the same dimensions but with a
        reduced number of  *maxpoints* points is returned.

        .. Note::

           Assumes that the first column is time.

        """
        return self._decimate(numkit.timeseries.mean_histogrammed_function, a, maxpoints, **kwargs)

    def decimate_circmean(self, a, maxpoints, **kwargs):
        """Return data *a* circmean-decimated on *maxpoints*.

        Histograms each column into *maxpoints* bins and calculates
        the weighted circular mean in each bin as the decimated data,
        using
        :func:`numkit.timeseries.circmean_histogrammed_function`. The
        coarse grained time in the first column contains the centers
        of the histogram time.

        If *a* contains <= *maxpoints* then *a* is simply returned;
        otherwise a new array of the same dimensions but with a
        reduced number of  *maxpoints* points is returned.

        Keywords *low* and *high* can be used to set the
        boundaries. By default they are [-pi, +pi].

        This method returns a **masked** array where jumps are flagged
        by an insertion of a masked point.

        .. Note::

           Assumes that the first column is time and that the data are
           in **degrees**.

        .. Warning::

           Breaking of arrays only works properly with a two-column
           array because breaks are only inserted in the x-column
           (a[0]) where y1 = a[1] has a break.

        """
        a_rad = numpy.vstack((a[0], numpy.deg2rad(a[1:])))
        b = self._decimate(numkit.timeseries.circmean_histogrammed_function, a_rad, maxpoints, **kwargs)
        y_ma, x_ma = break_array(b[1], threshold=numpy.pi, other=b[0])
        v = [y_ma]
        for y in b[2:]:
            v.append(break_array(y, threshold=numpy.pi)[0])
            if v[-1].shape != v[0].shape:
                raise ValueError("y dimensions have different breaks: you MUST deal with them separately")
        return numpy.vstack((x_ma, numpy.rad2deg(v)))

    def decimate_min(self, a, maxpoints, **kwargs):
        """Return data *a* min-decimated on *maxpoints*.

        Histograms each column into *maxpoints* bins and calculates
        the minimum in each bin as the decimated data, using
        :func:`numkit.timeseries.min_histogrammed_function`. The coarse grained
        time in the first column contains the centers of the histogram
        time.

        If *a* contains <= *maxpoints* then *a* is simply returned;
        otherwise a new array of the same dimensions but with a
        reduced number of  *maxpoints* points is returned.

        .. Note::

           Assumes that the first column is time.

        """
        return self._decimate(numkit.timeseries.min_histogrammed_function, a, maxpoints, **kwargs)

    def decimate_max(self, a, maxpoints, **kwargs):
        """Return data *a* max-decimated on *maxpoints*.

        Histograms each column into *maxpoints* bins and calculates
        the maximum in each bin as the decimated data, using
        :func:`numkit.timeseries.max_histogrammed_function`. The coarse grained
        time in the first column contains the centers of the histogram
        time.

        If *a* contains <= *maxpoints* then *a* is simply returned;
        otherwise a new array of the same dimensions but with a
        reduced number of  *maxpoints* points is returned.

        .. Note::

           Assumes that the first column is time.

        """
        return self._decimate(numkit.timeseries.max_histogrammed_function, a, maxpoints, **kwargs)

    def decimate_rms(self, a, maxpoints, **kwargs):
        """Return data *a* rms-decimated on *maxpoints*.

        Histograms each column into *maxpoints* bins and calculates
        the root mean square sum in each bin as the decimated data,
        using :func:`numkit.timeseries.rms_histogrammed_function`. The coarse
        grained time in the first column contains the centers of the
        histogram time.

        If *a* contains <= *maxpoints* then *a* is simply returned;
        otherwise a new array of the same dimensions but with a
        reduced number of  *maxpoints* points is returned.

        .. Note::

           Assumes that the first column is time.

        """
        return self._decimate(numkit.timeseries.rms_histogrammed_function, a, maxpoints, **kwargs)

    def decimate_percentile(self, a, maxpoints, **kwargs):
        """Return data *a* percentile-decimated on *maxpoints*.

        Histograms each column into *maxpoints* bins and calculates
        the percentile *per* in each bin as the decimated data, using
        :func:`numkit.timeseries.percentile_histogrammed_function`. The
        coarse grained time in the first column contains the centers
        of the histogram time.

        If *a* contains <= *maxpoints* then *a* is simply returned;
        otherwise a new array of the same dimensions but with a
        reduced number of  *maxpoints* points is returned.

        .. Note::

           Assumes that the first column is time.

        :Keywords:

        *per*
            percentile as a percentage, e.g. 75 is the value that splits
            the data into the lower 75% and upper 25%; 50 is the median [50.0]

        .. SeeAlso:: :func:`numkit.timeseries.regularized_function` with :func:`scipy.stats.scoreatpercentile`
        """
        return self._decimate(numkit.timeseries.percentile_histogrammed_function, a, maxpoints, **kwargs)

    def decimate_error(self, a, maxpoints, **kwargs):
        """Return data *a* error-decimated on *maxpoints*.

        Histograms each column into *maxpoints* bins and calculates an
        error estimate in each bin as the decimated data, using
        :func:`numkit.timeseries.error_histogrammed_function`. The
        coarse grained time in the first column contains the centers
        of the histogram time.

        If *a* contains <= *maxpoints* then *a* is simply returned;
        otherwise a new array of the same dimensions but with a
        reduced number of  *maxpoints* points is returned.

        .. SeeAlso:: :func:`numkit.timeseries.tcorrel`

        .. Note::

           Assumes that the first column is time.

           Does not work very well because often there are too few
           datapoints to compute a good autocorrelation function.

        """
        warnings.warn("Using undocumented decimate_error() is highly EXPERIMENTAL",
                      category=LowAccuracyWarning)
        return self._decimate(numkit.timeseries.error_histogrammed_function, a, maxpoints, **kwargs)

    def _decimate(self, func, a, maxpoints, **kwargs):
        ny = a.shape[-1]   # assume 2D array with last dimension varying fastest
        out = numpy.zeros((a.shape[0], maxpoints), dtype=float)

        t = a[0]
        for i in xrange(1, a.shape[0]):
            # compute regularised data for each column separately
            out[i], out[0] = func(t, a[i], bins=maxpoints, **kwargs)

        if maxpoints == self.maxpoints_default:  # only warn if user did not set maxpoints
            warnings.warn("Plot had %d datapoints > maxpoints = %d; decimated to %d regularly "
                          "spaced points from the histogrammed data with %s()."
                          % (ny, maxpoints, maxpoints, func.func_name),
                          category=AutoCorrectionWarning)
        return out


    def decimate_smooth(self, a, maxpoints, window="flat"):
        """Return smoothed data *a* decimated on approximately *maxpoints* points.

        1. Produces a smoothed graph using the smoothing window *window*;
           "flat" is a running average.
        2. select points at a step size approximatelt producing maxpoints

        If *a* contains <= *maxpoints* then *a* is simply returned;
        otherwise a new array of the same dimensions but with a
        reduced number of points (close to *maxpoints*) is returned.

        .. Note::

           Assumes that the first column is time (which will *never*
           be smoothed/averaged), except when the input array *a* is
           1D and therefore to be assumed to be data at equidistance
           timepoints.

        TODO:
        - Allow treating the 1st column as data
        """
        ny = a.shape[-1]   # assume 1D or 2D array with last dimension varying fastest
        # reduce size by averaging oover stepsize steps and then just
        # picking every stepsize data points.  (primitive --- can
        # leave out bits at the end or end up with almost twice of
        # maxpoints)
        stepsize = int(ny / maxpoints)
        if stepsize % 2 == 0:
            stepsize += 1  # must be odd for the running average/smoothing window
        out = numpy.empty_like(a)

        # smoothed
        out[0,:] = a[0]
        for i in xrange(1, a.shape[0]):
            # process columns because smooth() only handles 1D arrays :-p
            out[i,:] = numkit.timeseries.smooth(a[i], stepsize, window=window)

        if maxpoints == self.maxpoints_default:  # only warn if user did not set maxpoints
            warnings.warn("Plot had %d datapoints > maxpoints = %d; decimated to %d regularly "
                          "spaced points with smoothing (%r) over %d steps."
                          % (ny, maxpoints, ny/stepsize, window, stepsize),
                          category=AutoCorrectionWarning)
        return out[..., ::stepsize]

    def __getstate__(self):
        """custom pickling protocol: http://docs.python.org/library/pickle.html

        If :attr:`XVG.savedata` is ``False`` then any attributes in
        :attr:`XVG.__pickle_excluded` are *not* pickled as they are but simply
        pickled with the default value.
        """
        if self.savedata:
            d = self.__dict__
        else:
            # do not pickle the big array cache
            mangleprefix = '_'+self.__class__.__name__
            def demangle(k):
                """_XVG__array --> __array"""
                if k.startswith(mangleprefix):
                    k = k.replace(mangleprefix,'')
                return k
            d = {}
            for k in self.__dict__:
                d[k] = self.__pickle_excluded.get(demangle(k), self.__dict__[k])
        return d

    def __setstate__(self, d):
        # compatibility with older (pre 0.1.13) pickled instances
        if not 'savedata' in d:
            wmsg = "Reading pre 0.1.13 pickle file: setting savedata=False"
            warnings.warn(wmsg, category=DeprecationWarning)
            self.logger.warn(wmsg)
            d['savedata'] = False  # new default
        self.__dict__.update(d)


def break_array(a, threshold=numpy.pi, other=None):
    """Create a array which masks jumps >= threshold.

    Extra points are inserted between two subsequent values whose
    absolute difference differs by more than threshold (default is
    pi).

    Other can be a secondary array which is also masked according to
    *a*.

    Returns (*a_masked*, *other_masked*) (where *other_masked* can be
    ``None``)
    """
    assert len(a.shape) == 1, "Only 1D arrays supported"

    if other is not None:
        if not a.shape == other.shape:
            raise ValueError("arrays must be of identical shape")

    # jump occurs after the index in break
    breaks = numpy.where(numpy.abs(numpy.diff(a)) >= threshold)[0]
    # insert a blank after
    breaks += 1

    # is this needed?? -- no, but leave it here as a reminder
    #f2 = numpy.diff(a, 2)
    #up = (f2[breaks - 1] >= 0)  # >0: up, <0: down
    # sort into up and down breaks:
    #breaks_up = breaks[up]
    #breaks_down = breaks[~up]

    # new array b including insertions for all the breaks
    m = len(breaks)
    b = numpy.empty((len(a) + m))
    # calculate new indices for breaks in b, taking previous insertions into account
    b_breaks = breaks + numpy.arange(m)
    mask =  numpy.zeros_like(b, dtype=numpy.bool)
    mask[b_breaks] = True
    b[~mask] = a
    b[mask] = numpy.NAN

    if other is not None:
        c = numpy.empty_like(b)
        c[~mask] = other
        c[mask] = numpy.NAN
        ma_c = numpy.ma.array(c, mask=mask)
    else:
        ma_c = None

    return numpy.ma.array(b, mask=mask), ma_c




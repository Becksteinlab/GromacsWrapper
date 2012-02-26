# GromacsWrapper: formats.py
# Copyright (c) 2009-2011 Oliver Beckstein <orbeckst@gmail.com>
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

.. autoclass:: XVG
   :members:
"""


from __future__ import with_statement
import os, errno
import re
import warnings

import numpy

from gromacs import ParseError, AutoCorrectionWarning
import gromacs.utilities as utilities
from gromacs.odict import odict

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

    .. Note:: - Only simple XY or NXY files are currently supported, *not*
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

    def __init__(self, filename=None, names=None, permissive=False, **kwargs):
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

        """
        self.__array = None          # cache for array (BIG) (used by XVG.array)
        self.__cache = {}            # cache for computed results
        self.savedata = kwargs.pop('savedata', False)
        if not filename is None:
            self._init_filename(filename)  # note: reading data from file is delayed until required
        if names is None:
            self.names = []
        else:
            try:
                self.names = names.split(',')
            except AttributeError:
                self.names = names
        self.permissive = permissive
        self.stride = kwargs.pop('stride', 1)
        self.corrupted_lineno = None      # must parse() first before this makes sense
        # default number of data points for calculating correlation times via FFT
        self.ncorrel = kwargs.pop('ncorrel', 25000)
        self.__correlkwargs = {}          # see set_correlparameters()

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
                if line.startswith(('#', '@')) or len(line) == 0:
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
          *kwargs*
               All other keyword arguments are passed on to :func:`pylab.plot`.
        """
        import pylab

        columns = kwargs.pop('columns', Ellipsis)         # slice for everything
        maxpoints = kwargs.pop('maxpoints', self.maxpoints_default)
        transform = kwargs.pop('transform', lambda x: x)  # default is identity transformation
        method = kwargs.pop('method', "mean")

        # (decimate/smooth o slice o transform)(array)
        a = self.decimate(method, numpy.asarray(transform(self.array))[columns], maxpoints)

        if len(a.shape) == 1:
            # special case: plot against index; plot would do this automatically but
            # we'll just produce our own xdata and pretend that this was X all along
            X = numpy.arange(len(a))
            a = numpy.concatenate([[X], [a]])  # does NOT overwrite original a but make a new one

        # now deal with infs, nans etc AFTER all transformations (needed for plotting across inf/nan)
        ma = numpy.ma.MaskedArray(a, mask=numpy.logical_not(numpy.isfinite(a)))

        # finally plot
        kwargs['xdata'] = ma[0]          # abscissa set separately
        pylab.plot(ma[1:].T, **kwargs)   # plot all other columns in parallel

    def errorbar(self, **kwargs):
        """Quick hack: errorbar plot.

        Set *columns* keyword to select [x, y, dy] or [x, y, dx, dy],
        e.g. ``columns=[0,1,2]``. See :meth:`XVG.plot` for details.
        """
        # TODO: This was copy&paste+modify from plot() and hence most of the code is duplicated
        #       -- needs to be done properly!!

        import pylab

        kwargs.setdefault('capsize', 0)
        kwargs.setdefault('elinewidth', 1)
        kwargs.setdefault('alpha', 0.3)
        kwargs.setdefault('fmt', None)

        columns = kwargs.pop('columns', Ellipsis)         # slice for everything
        maxpoints = kwargs.pop('maxpoints', self.maxpoints_default)
        transform = kwargs.pop('transform', lambda x: x)  # default is identity transformation
        method = kwargs.pop('method', "mean")

        # (decimate/smooth o slice o transform)(array)
        a = self.decimate(method, numpy.asarray(transform(self.array))[columns], maxpoints)

        if len(a.shape) == 1:
            # special case: plot against index; plot would do this automatically but
            # we'll just produce our own xdata and pretend that this was X all along
            X = numpy.arange(len(a))
            a = numpy.concatenate([[X], [a]])  # does NOT overwrite original a but make a new one

        # now deal with infs, nans etc AFTER all transformations (needed for plotting across inf/nan)
        ma = numpy.ma.MaskedArray(a, mask=numpy.logical_not(numpy.isfinite(a)))

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

        pylab.errorbar(X, Y, **kwargs)

    def decimate(self, method, a, maxpoints, **kwargs):
        """Decimate data *a* to *maxpoints* using *method*.

        Methods:

          * "mean", uses :meth:`XVG.decimate_mean` to coarse grain by
            averaging the data in bins along the time axis

          * "smooth", uses :meth:`XVG.decimate_smooth` to subsample
            from a smoothed function (generated with a running average
            of the coarse graining step size derived from the original
            number of data points and *maxpoints*)
        """
        methods = {'mean': self.decimate_mean,
                   'smooth': self.decimate_smooth,
                   }
        return methods[method](a, maxpoints, **kwargs)

    def decimate_mean(self, a, maxpoints, **kwargs):
        """Return data *a* decimated on *maxpoints*.

        Histograms each column into *maxpoints* bins and calculates
        the weighted average in each bin as the decimated data, using
        :func:`numkit.timeseries.mean_histogram`. The coarse grained
        time in the first column contains the centers of the histogram
        time.

        If *a* contains <= *maxpoints* then *a* is simply returned;
        otherwise a new array of the same dimensions but with a
        reduced number of  *maxpoints* points is returned.

        .. Warning::

           Assumes that the first column is time, except when the
           input array *a* is 1D and therefore to be assumed to be
           data at equidistance timepoints.

        """
        if len(a.shape) == 1:
            # add first column as index
            # (probably should do this in class/init anyway...)
            X = numpy.arange(len(a))
            a = numpy.vstack([X, a])
        ny = a.shape[-1]   # assume 1D or 2D array with last dimension varying fastest
        if maxpoints is None or ny <= maxpoints:
            return a

        out = numpy.zeros((a.shape[0], maxpoints), dtype=float)

        t = a[0]
        for i in xrange(1, a.shape[0]):
            # compute regularised data for each column separately
            out[i], out[0] = numkit.timeseries.mean_histogrammed_function(t, a[i], bins=maxpoints)

        if maxpoints == self.maxpoints_default:  # only warn if user did not set maxpoints
            warnings.warn("Plot had %d datapoints > maxpoints = %d; decimated to %d regularly "
                          "spaced points by computing the bin-means from the histogrammed data."
                          % (ny, maxpoints, maxpoints),
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

        .. Warning::

           Assumes that the first column is time (which will *never*
           be smoothed/averaged), except when the input array *a* is
           1D and therefore to be assumed to be data at equidistance
           timepoints.

        TODO:
        - Allow treating the 1st column as data
        """

        ny = a.shape[-1]   # assume 1D or 2D array with last dimension varying fastest
        if maxpoints is None or ny <= maxpoints:
            return a

        # reduce size by averaging oover stepsize steps and then just
        # picking every stepsize data points.  (primitive --- can
        # leave out bits at the end or end up with almost twice of
        # maxpoints)
        stepsize = int(ny / maxpoints)
        if stepsize % 2 == 0:
            stepsize += 1  # must be odd for the running average/smoothing window
        out = numpy.empty_like(a)

        # smoothed
        if len(a.shape) == 1:
            out[:] = numkit.timeseries.smooth(a, stepsize, window=window)
        else:
            out[0,:] = a[0]
            for i in xrange(1, a.shape[0]):
                # process columns because smooth() only handles 1D arrays :-p
                # should change y in place
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



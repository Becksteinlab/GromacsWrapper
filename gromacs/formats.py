# GromacsWrapper: formats.py
# Copyright (c) 2009-2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
:mod:`gromacs.formats` -- Accessing various files
=================================================

This module contains classes that represent data files on
disk. Typically one creates an instance and

- reads from a file using a :meth:`read` method, or

- populates the instance (in the simplest case with a :meth:`set`
  method) and the uses the :meth:`write` method to write the data to
  disk in the appropriate format.

For function data there typically also exists a :meth:`plot` method
which produces a graph (using matplotlib).

The module defines some classes that are used in other modules; they
do *not* make use of :mod:`gromacs.tools` or :mod:`gromacs.cbook` and
can be safely imported at any time.


Classes
-------

.. autoclass:: XVG
   :members:
.. autoclass:: NDX
   :members:
.. autoclass:: uniqueNDX
   :members:
.. autoclass:: GRO
   :members:

   (Not implemented yet)
"""
from __future__ import with_statement

__docformat__ = "restructuredtext en"

import os
import re
import warnings
import errno
import operator

import numpy

from odict import odict

import utilities
from gromacs import ParseError, AutoCorrectionWarning

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
    logger = logging.getLogger('gromacs.formats.XVG')  # for pickling: must be class-level

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
        """
        from numkit.timeseries import tcorrel
        from gromacs.analysis.collections import Collection
        t = self.array[0,::nstep]
        r = Collection([tcorrel(t, Y, nstep=1, **kwargs) for Y in self.array[1:,::nstep]])
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

        .. _p526: http://books.google.co.uk/books?id=XmyO2oRUg0cC&lpg=PP1&ots=Zw5FY2j2DS&dq=Frenkel%20and%20Smit%2C&pg=PA526
        """
        return self._correlprop('sigma')

    @property
    def tc(self):
        """Correlation time of the data.

        See :meth:`XVG.error` for details.
        """
        return self._correlprop('tc')

    def parse(self):
        """Read and cache the file as a numpy array.

        The array is returned with column-first indexing, i.e. for a data file with
        columns X Y1 Y2 Y3 ... the array a will be a[0] = X, a[1] = Y1, ... .
        """
        self.corrupted_lineno = []
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
                ncol = len(row)
                rows.append(row)
        try:
            self.__array = numpy.array(rows).transpose()    # cache result
        except:
            self.logger.error("%s: Failed reading XVG file, possibly data corrupted. "
                              "Check the last line of the file...", self.real_filename)
            raise

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
               subsample the data at regular intervals. [10000]
          *kwargs*
               All other keyword arguments are passed on to :func:`pylab.plot`.
        """
        import pylab

        maxpoints_default = 10000
        columns = kwargs.pop('columns', Ellipsis)         # slice for everything
        maxpoints = kwargs.pop('maxpoints', maxpoints_default)
        transform = kwargs.pop('transform', lambda x: x)  # default is identity transformation
        a = numpy.asarray(transform(self.array))[columns] # (slice o transform)(array)

        ny = a.shape[-1]   # assume 1D or 2D array with last dimension varying fastest
        if not maxpoints is None and ny > maxpoints:
            # reduce size by subsampling (primitive --- can leave out
            # bits at the end or end up with almost twice of maxpoints)
            stepsize = int(ny / maxpoints)
            a = a[..., ::stepsize]
            if maxpoints == maxpoints_default:  # only warn if user did not set maxpoints
                warnings.warn("Plot had %d datapoints > maxpoints = %d; subsampled to %d regularly spaced points." 
                              % (ny, maxpoints, a.shape[-1]), category=AutoCorrectionWarning)

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
        
        Set columns to select [x, y, dy].
        """
        import pylab

        kwargs.setdefault('capsize', 0)
        kwargs.setdefault('elinewidth', 1)
        kwargs.setdefault('alpha', 0.3)
        kwargs.setdefault('fmt', None)

        maxpoints_default = 10000
        columns = kwargs.pop('columns', Ellipsis)         # slice for everything
        maxpoints = kwargs.pop('maxpoints', maxpoints_default)
        transform = kwargs.pop('transform', lambda x: x)  # default is identity transformation
        a = numpy.asarray(transform(self.array))[columns] # (slice o transform)(array)

        ny = a.shape[-1]   # assume 1D or 2D array with last dimension varying fastest
        if not maxpoints is None and ny > maxpoints:
            # reduce size by subsampling (primitive --- can leave out
            # bits at the end or end up with almost twice of maxpoints)
            stepsize = int(ny / maxpoints)
            a = a[..., ::stepsize]
            if maxpoints == maxpoints_default:  # only warn if user did not set maxpoints
                warnings.warn("Plot had %d datapoints > maxpoints = %d; subsampled to %d regularly spaced points." 
                              % (ny, maxpoints, a.shape[-1]), category=AutoCorrectionWarning)

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
            kwargs['yerr'] = ma[2]

        pylab.errorbar(X, Y, **kwargs)

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


class NDX(odict, utilities.FileUtils):
    """Gromacs index file.

    Represented as a ordered dict where the keys are index group names and
    values are numpy arrays of atom numbers.

    Use the :meth:`NDX.read` and :meth:`NDX.write` methods for
    I/O. Access groups by name via the :meth:`NDX.get` and
    :meth:`NDX.set` methods.

    Alternatively, simply treat the :class:`NDX` instance as a
    dictionary. Setting a key automatically transforms the new value
    into a integer 1D numpy array (*not* a set, as would be the
    :program:`make_ndx` behaviour).

    .. Note:: The index entries themselves are ordered and can contain 
              duplicates so that output from NDX can be easily used for 
              :program:`g_dih` and friends. If you need set-like behaviour
              you will have do use :class:`gromacs.formats.uniqueNDX` or
              :class:`gromacs.cbook.IndexBuilder` (which uses
              :program:`make_ndx` throughout).

    **Example**

      Read index file, make new group and write to disk::

        ndx = NDX()
        ndx.read('system.ndx')
        print ndx['Protein']       
        ndx['my_group'] = [2, 4, 1, 5]   # add new group
        ndx.write('new.ndx')

      Or quicker (replacing the input file ``system.ndx``)::

        ndx = NDX('system')          # suffix .ndx is automatically added
        ndx['chi1'] = [2, 7, 8, 10]
        ndx.write()

    """
    default_extension = "ndx"
    
    # match:  [ index_groupname ]
    SECTION = re.compile("""\s*\[\s*(?P<name>\S.*\S)\s*\]\s*""")

    #: standard ndx file format: 15 columns
    ncol = 15
    #: standard ndx file format: '%6d'
    format = '%6d'

    def __init__(self, filename=None, **kwargs):
        super(NDX, self).__init__(**kwargs)  # can use kwargs to set dict! (but no sanity checks!)

        if not filename is None:
            self._init_filename(filename)
            self.read(filename)

    def read(self, filename=None):
        """Read and parse index file *filename*."""        
        self._init_filename(filename)
        
        data = odict()
        with open(self.real_filename) as ndx:
            current_section = None
            for line in ndx:
                line = line.strip()
                if len(line) == 0:
                    continue
                m = self.SECTION.match(line)
                if m:
                    current_section = m.group('name')
                    data[current_section] = []  # can fail if name not legal python key
                    continue
                if not current_section is None:
                    data[current_section].extend(map(int, line.split()))

        super(NDX,self).update(odict([(name, self._transform(atomnumbers))
                                     for name, atomnumbers in data.items()]))

    def write(self, filename=None, ncol=ncol, format=format):
        """Write index file to *filename* (or overwrite the file that the index was read from)"""
        with open(self.filename(filename, ext='ndx'), 'w') as ndx:
            for name in self:
                atomnumbers = self._getarray(name)  # allows overriding
                ndx.write('[ %s ]\n' % name)
                for k in xrange(0, len(atomnumbers), ncol):
                    line = atomnumbers[k:k+ncol].astype(int)   # nice formatting in ncol-blocks
                    n = len(line)
                    ndx.write((" ".join(n*[format])+'\n') % tuple(line))
                ndx.write('\n')

    def get(self, name):
        """Return index array for index group *name*."""
        return self[name]

    def set(self, name, value):
        """Set or add group *name* as a 1D numpy array."""
        self[name] = value

    def size(self, name):
        """Return number of entries for group *name*."""
        return len(self[name])

    @property
    def groups(self):
        """Return a list of all groups."""
        return self.keys()

    @property
    def sizes(self):
        """Return a dict with group names and number of entries,"""
        return dict([(name, len(atomnumbers)) for name, atomnumbers in self.items()])

    @property
    def ndxlist(self):
        """Return a list of groups in the same format as  :func:`gromacs.cbook.get_ndx_groups`.

        Format:
           [ {'name': group_name, 'natoms': number_atoms, 'nr':  # group_number}, ....]
        """
        return [{'name': name, 'natoms': len(atomnumbers), 'nr': nr+1} for
                nr,(name,atomnumbers) in enumerate(self.items())]

    def _getarray(self, name):
        """Helper getter that is used in write(). 
        Override when using a _transform that stores something that
        cannot be indexed, e.g. when using set()s.
        """
        return self[name]

    def _transform(self, v):
        """Transform input to the stored representation.
        
        Override eg with ``return set(v)`` for index lists as sets.
        """
        return numpy.ravel(v).astype(int)
        
    def __setitem__(self, k, v):
        super(NDX, self).__setitem__(k, self._transform(v))

    def setdefault(*args,**kwargs):
        raise NotImplementedError


class IndexSet(set):
    """set which defines '+' as union (OR) and '-' as intersection  (AND)."""
    def __add__(self, x):
        return self.union(x)
    def __sub__(self, x):
        return self.intersection(x)


class uniqueNDX(NDX):
    """Index that behaves like make_ndx, i.e. entries behaves as sets,
    not lists.

    The index lists behave like sets:
    - adding sets with '+' is equivalent to a logical OR: x + y == "x | y"
    - subtraction '-' is AND: x - y == "x & y"
    - see :meth:`~gromacs.formats.join` for ORing multiple groups (x+y+z+...)

    **Example** ::
       I = uniqueNDX('system.ndx')
       I['SOLVENT'] = I['SOL'] + I['NA+'] + I['CL-']    
    """

    def join(self, *groupnames):
        """Return an index group that contains atoms from all  *groupnames*.

        The method will silently ignore any groups that are not in the
        index.

        **Example**

        Always make a solvent group from water and ions, even if not
        all ions are present in all simulations::
        
           I['SOLVENT'] = I.join('SOL', 'NA+', 'K+', 'CL-')        
        """
        return self._sum([self[k] for k in groupnames if k in self])

    def _sum(self, sequence):
        return reduce(operator.add, sequence)

    def _transform(self, v):
        return IndexSet(v)

    def _getarray(self, k):
        return numpy.sort(numpy.fromiter(self[k],dtype=int,count=len(self[k])))

    

# or use list of these?
# class IndexGroup(dict):
#     def __init__(self, groupnumber=None, name="empty", atomnumbers=None, **kwargs):
#         atomnumbers = atomnumbers or []
#         _atomnumbers = numpy.asarray(atomnumbers).astype(int)
#         super(IndexGroup, self).__init__(name=str(name),
#                                          atomnumbers=_atomnumbers,
#                                          nr=groupnumber)
    
class GRO(utilities.FileUtils):
    """Class that represents a GROMOS (gro) structure file.


    File format:
    """
    default_extension = "gro"
    logger = logging.getLogger('gromacs.formats.GRO')

    def __init__(self, **kwargs):

        raise NotImplementedError
        
        filename = kwargs.pop('filename',None)
        super(GRO, self).__init__(**kwargs)

        if not filename is None:
            self._init_filename(filename)
            self.read(filename)

    def read(self, filename=None):
        """Read and parse index file *filename*."""        
        self._init_filename(filename)
        
        with open(self.real_filename) as gro:
            pass
        


class MDP(odict, utilities.FileUtils):
    """Class that represents a Gromacs mdp run input file.

    The MDP instance is an ordered dictionary.

      - *Parameter names* are keys in the dictionary.
      - *Comments* are sequentially numbered with keys Comment0001,
        Comment0002, ...
      - *Empty lines* are similarly preserved as Blank0001, ....

    When writing, the dictionary is dumped in the recorded order to a
    file. Inserting keys at a specific position is not possible.

    Currently, comments after a parameter on the same line are
    discarded. Leading and trailing spaces are always stripped.

    .. SeeAlso:: For editing a mdp file one can also use
                :func:`gromacs.cbook.edit_mdp` (which works like a
                poor replacement for sed).
    """
    default_extension = "mdp"
    logger = logging.getLogger('gromacs.formats.MDP')

    COMMENT = re.compile("""\s*;\s*(?P<value>.*)""")   # eat initial ws
    # see regex in cbook.edit_mdp()
    PARAMETER = re.compile("""
                            \s*(?P<parameter>[^=]+?)\s*=\s*  # parameter (ws-stripped), before '='
                            (?P<value>[^;]*)                # value (stop before comment=;)
                            (?P<comment>\s*;.*)?            # optional comment           
                            """, re.VERBOSE)
        
    def __init__(self, filename=None, autoconvert=True, **kwargs):
        """Initialize mdp structure.

        :Arguments:
          *filename*
              read from mdp file
          *autoconvert* : boolean
              ``True`` converts numerical values to python numerical types;
              ``False`` keeps everything as strings [``True``]
          *kwargs*
              Populate the MDP with key=value pairs. (NO SANITY CHECKS; and also
              does not work for keys that are not legal python variable names such
              as anything that includes a minus '-' sign or starts with a number).
        """
        super(MDP, self).__init__(**kwargs)  # can use kwargs to set dict! (but no sanity checks!)

        self.autoconvert = autoconvert
        
        if not filename is None:
            self._init_filename(filename)
            self.read(filename)

    def _transform(self, value):
        if self.autoconvert:
            return autoconvert(value)
        else:
            return value

    def read(self, filename=None):
        """Read and parse mdp file *filename*."""        
        self._init_filename(filename)

        def BLANK(i):
            return "B%04d" % i
        def COMMENT(i):
            return "C%04d" % i
        
        data = odict()
        iblank = icomment = 0
        with open(self.real_filename) as mdp:
            for line in mdp:
                line = line.strip()
                if len(line) == 0:
                    iblank += 1
                    data[BLANK(iblank)] = ''
                    continue
                m = self.COMMENT.match(line)
                if m:
                    icomment += 1
                    data[COMMENT(icomment)] = m.group('value')
                    continue
                # parameter
                m = self.PARAMETER.match(line)
                if m:
                    # check for comments after parameter?? -- currently discarded
                    parameter = m.group('parameter')
                    value =  self._transform(m.group('value'))
                    data[parameter] = value
                else:
                    errmsg = '%(filename)r: unknown line in mdp file, %(line)r' % vars()
                    self.logger.error(errmsg)
                    raise ParseError(errmsg)

        super(MDP,self).update(data)


    def write(self, filename=None, skipempty=False):
        """Write mdp file to *filename*.

        :Keywords:
           *filename*
               output mdp file; default is the filename the mdp
               was read from
           *skipempty* : boolean
               ``True`` removes any parameter lines from output that
               contain empty values [``False``]

        .. Note:: Overwrites the file that the mdp was read from if no
                  *filename* supplied.
        """

        with open(self.filename(filename, ext='mdp'), 'w') as mdp:
            for k,v in self.items():
                if k[0] == 'B':        # blank line
                    mdp.write("\n")
                elif k[0] == 'C':      # comment
                    mdp.write("; %(v)s\n" % vars())
                else:                  # parameter = value
                    if skipempty and (v == '' or v is None):
                        continue
                    mdp.write("%(k)s = %(v)s\n" % vars())


def autoconvert(s):
    """Convert input to a numerical type if possible.

    1. A non-string object is returned as it is
    2. Try conversion to int, float, str.
    """
    if not type(s) is str:
        return s
    for converter in int, float, str:   # try them in increasing order of lenience
        try:
            return converter(s)
        except ValueError:
            pass
    raise ValueError("Failed to autoconvert %r" % s)
                

# $Id$
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
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


"""
from __future__ import with_statement

__docformat__ = "restructuredtext en"

import os
import re
import warnings
import errno

import numpy

import utilities

class XVG(utilities.FileUtils):
    """Class that represents the numerical data in a grace xvg file.

    All data must be numerical. :const:`NAN` and :const:`INF` values are
    supported via python's :func:`float` builtin function.

    The :attr:`~XVG.array` attribute can be used to access the
    the array once it has been read and parsed.

    Conceptually, the file on disk and the XVG instance are considered the same
    data. This means that whenever the filename for I/O (:meth:`XVG.read` and
    :meth:`XVG.write`) is changed then the filename associated with the
    instance is also changed to reflect the association between file and
    instance.

    .. Note:: - Only simple XY or NXY files are currently supported, not
                Grace files that contain multiple data sets separated by '&'.
              - Any kind of formatting (xmgrace commands) are discarded.
    """
    def __init__(self, filename=None):
        """Initialize the class from a xvg file.

        :Arguments: *filename* is the xvg file; it can only be of type XY or
                    NXY. If it is supplied then it is read and parsed when
                    :attr:`XVG.array` is accessed.
        """
        self.__array = None          # cache for array property
        if not filename is None:
            self._init_filename(filename)  # reading is delayed until required

    def _init_filename(self, filename):
        """Initialize the current filename :attr:`XVG.real_filename` of the object.

        - The first invocation must have ``filename != None``; this will set a
          default filename with suffix ".xvg" unless another one was supplied.
        - Subsequent invocations either change the filename accordingly or
          ensure that the default filename is set with the proper suffix.        
        """
        f = self.filename(filename, ext='xvg', use_my_ext=True, set_default=True)
        #: Current full path of the object for reading and writing I/O.
        self.real_filename = os.path.realpath(f)

    def read(self, filename=None):
        """Read and parse xvg file *filename*."""
        self._init_filename(filename)
        self.parse()

    def write(self, filename=None):
        """Write array to xvg file *filename* in NXY format."""
        self._init_filename(filename)
        with open(self.real_filename, 'w') as xvg:
            xvg.write("# xmgrace compatible NXY data file\n"
                      "# Written by gromacs.formats.XVG()\n")
            for xyy in self.array.T:
                xyy.tofile(xvg, sep=" ", format="%-8s")     # quick and dirty ascii output...
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
        
    def parse(self):
        """Read and cache the file as a numpy array.

        The array is returned with column-first indexing, i.e. for a data file with
        columns X Y1 Y2 Y3 ... the array a will be a[0] = X, a[1] = Y1, ... .
        """
        # cannot use numpy.loadtxt() because xvg can have two types of 'comment' lines
        with open(self.real_filename) as xvg:
            rows = []
            for line in xvg:
                line = line.strip()
                if line.startswith(('#', '@')) or len(line) == 0:
                    continue
                if line.startswith('&'):
                    raise NotImplementedError('Sorry only simple NXY format is supported.')
                rows.append(map(float, line.split()))
        self.__array = numpy.array(rows).transpose()    # cache result

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
            if maxpoints == maxpoints_default:
                warnings.warn("Plot had %d datapoints > maxpoints = %d; subsampled to %d regularly spaced points." 
                              % (ny, maxpoints, a.shape[-1]), category=AutoCorrectionWarning)

        if len(a.shape) == 1:
            # special case: plot against index; plot would do this automatically but 
            # we'll just produce our own xdata and pretend that this was X all along
            X = numpy.arange(len(a))
            a = numpy.concatenate([[X], [a]])  # does NOT overwrite original a but make a new one
        kwargs['xdata'] = a[0]          # abscissa set separately
        pylab.plot(a[1:].T, **kwargs)   # plot all other columns in parallel
        

class NDX(dict, utilities.FileUtils):
    """Gromacs index file.

    Represented as a dict where the keys are index group names and
    values are numpy arrays of atom numbers.

    Use the :meth:`NDX.read` and :meth:`NDX.write` methods for
    I/O. Access groups by name via the :meth:`NDX.get` and
    :meth:`NDX.set` methods.

    Alternatively, simply treat the NDX instance as a
    dictionary. Setting a key automatically transforms the new value
    into a integer 1D numpy array.

    Example:

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

    .. Note:: When writing an index file then the order in which groups are
              written to the filke is undetermined. This is arguably a bug and
              should be fixed in future releases.
    """
    # TODO: use a ordered dict to preserve order of groups (important when
    #       accessing groups by number)

    # match:  [ index_groupname ]
    SECTION = re.compile("""\s*\[\s*(?P<name>\S.*\S)\s*\]\s*""")

    #: standard ndx file format: 15 columns
    ncol = 15
    #: standard ndx file format: '%6d'
    format = '%6d'

    def __init__(self, **kwargs):
        super(NDX, self).__init__()

        filename = kwargs.pop('filename',None)
        if not filename is None:
            self._init_filename(filename)
            self.read()

    def _init_filename(self, filename=None):
        filename = self.filename(filename, ext='ndx')
        self.real_filename = os.path.realpath(filename)  # use full path for accessing data        

    def read(self, filename=None):
        """Read and parse index file *filename*."""        
        self._init_filename(filename)
        
        with open(self.real_filename) as ndx:
            data = {}
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

        super(NDX,self).update(dict([(name, numpy.array(atomnumbers))
                                     for name, atomnumbers in data.items()]))

    def write(self, filename=None, ncol=ncol, format=format):
        """Write index file to *filename* (or overwrite the file that the index was read from)"""
        with open(self.filename(filename, ext='ndx'), 'w') as ndx:
            for name, atomnumbers in self.items():
                ndx.write('[ %s ]\n' % name)
                for k in xrange(0, len(atomnumbers)/ncol + 1):
                    line = atomnumbers[k:k+ncol].astype(int)   # nice formatting in ncol-blocks
                    n = len(line)
                    ndx.write((" ".join(n*[format])+'\n') % tuple(line))
                ndx.write('\n')

    def get(self, name):
        """Return index array for idnex group *name*."""
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

    def __setitem__(self, k, v):
        super(NDX, self).__setitem__(k, numpy.ravel(v).astype(int))

    def update(self,*args,**kwargs):
        raise NotImplementedError

    def setdefault(*args,**kwargs):
        raise NotImplementedError
    
    
    
    

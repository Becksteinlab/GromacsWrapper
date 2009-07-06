# $Id$
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
``gromacs.utilities`` -- Helper functions and classes
=====================================================

The module defines some convenience functions and classes that are
used in other modules; they do *not* make use of :mod:`gromacs.tools`
or :mod:`gromacs.cbook` and can be safely imported.

Classes
-------

:class:`FileUtils` provides functions related to filename handling. It
can be used as a mixin class. The :class:`gromacs.analysis.Simulation`
class is derived from it.

.. autoclass:: FileUtils
   :members:
.. autoclass:: AttributeDict
.. autoclass:; XVG

Functions
---------

Some additiona convenience functions:

.. autofunction:: anyopen
.. autofunction:: iterable
.. autofunction:: convert_aa_code
.. autofunction:: unlink_f
.. autofunction:: unlink_gmx
.. autofunction:: unlink_gmx_backups

Data
----

.. autodata:: amino_acid_codes

"""
from __future__ import with_statement

__docformat__ = "restructuredtext en"

import os
import glob
import warnings
import errno
import bz2, gzip
import numpy

def Property(func):
    """Simple decorator wrapper to make full fledged properties.
    See eg http://adam.gomaa.us/blog/2008/aug/11/the-python-property-builtin/
    """
    return property(**func())
    

class AttributeDict(dict):
    """A dictionary with pythonic access to keys as attributes --- useful for interactive work."""
    def __getattribute__(self,x):
        try:
            return super(AttributeDict,self).__getattribute__(x)
        except AttributeError:
            return self[x]
    def __setattr__(self,name,value):
        try:
            super(AttributeDict,self).__setitem__(name, value)
        except KeyError:
            super(AttributeDict,self).__setattr__(name, value)            

def anyopen(datasource):
    """Open datasource and return a stream."""
    if hasattr(datasource,'next') or hasattr(datasource,'readline'):
        stream = datasource
        filename = '(%s)' % stream.name  # maybe that does not always work?        
    else:
        stream = None
        filename = datasource
        for openfunc in bz2.BZ2File, gzip.open, file:   # file should be last
            stream = _get_stream(datasource, openfunc)
            if not stream is None:
                break
        if stream is None:
            raise IOError("Cannot open %r for reading." % filename)
    return stream, filename

def _get_stream(filename, openfunction=file):
    try:
        stream = openfunction(filename,'r')
    except IOError:
        return None

    try:
        stream.readline()
        stream.close()
        stream = openfunction(filename,'r')
    except IOError:
        stream.close()
        stream = None
    return stream

#: translation table for 1-letter codes --> 3-letter codes 
amino_acid_codes = {'A':'ALA', 'C':'CYS', 'D':'ASP', 'E':'GLU',
                    'F':'PHE', 'G':'GLY', 'H':'HIS', 'I':'ILE',
                    'K':'LYS', 'L':'LEU', 'M':'MET', 'N':'ASN',
                    'P':'PRO', 'Q':'GLN', 'R':'ARG', 'S':'SER',
                    'T':'THR', 'V':'VAL', 'W':'TRP', 'Y':'TYR'}
inverse_aa_codes = dict([(three, one) for one,three in amino_acid_codes.items()])

def convert_aa_code(x):
    """Converts between 3-letter and 1-letter amino acid codes."""
    if len(x) == 1:
        return amino_acid_codes[x.upper()]
    elif len(x) == 3:
        return inverse_aa_codes[x.upper()]
    else:
        raise ValueError("Can only convert 1-letter or 3-letter amino acid codes, "
                         "not %r" % x)


class FileUtils(object):
    """Mixin class to provide additional IO capabilities."""

    def filename(self,filename=None,ext=None,set_default=False,use_my_ext=False):
        """Supply a file name for the class object.

        Typical uses::

           fn = filename()             ---> <default_filename>
           fn = filename('name.ext')   ---> 'name'
           fn = filename(ext='pickle') ---> <default_filename>'.pickle'
           fn = filename('name.inp','pdf') --> 'name.pdf'
           fn = filename('foo.pdf',ext='png',use_my_ext=True) --> 'foo.pdf'

        The returned filename is stripped of the extension
        (``use_my_ext=False``) and if provided, another extension is
        appended. Chooses a default if no filename is given.  

        Raises a ``ValueError`` exception if no default file name is known.

        If ``set_default=True`` then the default filename is also set.

        ``use_my_ext=True`` lets the suffix of a provided filename take
        priority over a default ``ext`` tension.
        """
        if filename is None:
            if not hasattr(self,'_filename'):
                self._filename = None        # add attribute to class 
            if self._filename:
                filename = self._filename
            else:
                raise ValueError("A file name is required because no default file name was defined.")
            my_ext = None
        else:
            filename, my_ext = os.path.splitext(filename)
            if set_default:                  # replaces existing default file name
                self._filename = filename
        if my_ext and use_my_ext:  
            ext = my_ext
        if ext is not None:
            if ext.startswith('.'):
                ext = ext[1:]  # strip a dot to avoid annoying mistakes
            filename = filename + '.' + ext
        return filename

    def check_file_exists(self, filename, resolve='exception'):
        """If a file exists then continue with the action specified in ``resolve``.

        ``resolve`` must be one of

        "ignore"
              always return ``False``
        "indicate"
              return ``True`` if it exists
        "warn"
              indicate and issue a :exc:`UserWarning`
        "exception"
              raise :exc:`IOError` if it exists
        """
        def _warn(x):
            warnings.warn("File %r already exists." % x)
            return True
        def _raise(x):
            raise IOError("File %r already exists." % x)
        solutions = {'ignore': lambda x: False,      # file exists, but we pretend that it doesn't
                     'indicate': lambda x: True,     # yes, file exists
                     'warn': _warn,
                     'warning': _warn,
                     'exception': _raise,
                     'raise': _raise,
                     }
        if not os.path.isfile(filename):
            return False
        else:
            return solutions[resolve](filename)


class XVG(FileUtils):
    """Class that represents a grace xvg file. Read-only at the moment."""
    def __init__(self, filename):
        self.filename = filename
        self.real_filename = os.path.realpath(filename)  # use full path for accessing data
        self.__array = None          # cache for array property

    def asarray(self):
        """Return data of the file as numpy array.

        The array is returned with column-first indexing, i.e. for a data file with
        columns X Y1 Y2 Y3 ... the array a will be a[0] = X, a[1] = Y1, ... .

        All data must be numerical. ``NAN`` and ``INF`` values are
        supported via python's :func:`float` builtin function.

        Instead of using this function one can also use the attr:`~XVG.array`
        attribute to access a cached version of the array.

        .. Note:: Only simple XY or NXY files are currently supported, not
                  Grace files that contain multiple data sets separated by '&'.
        """
        with open(self.real_filename) as xvg:
            rows = []
            for line in xvg:
                line = line.strip()
                if line.startswith(('#', '@')) or len(line) == 0:
                    continue
                if line.startswith('&'):
                    raise NotImplementedError('Sorry only simple NXY format is supported.')
                rows.append(map(float, line.split()))
        return numpy.array(rows).transpose()

    @Property
    def array():
        doc = """Represent xvg data as a (cached) numpy array. 
              See meth:`~XVG.asarray` for details."""
        def fget(self):
            if self.__array is None:
                self.__array = self.asarray()
            return self.__array
        return locals()

    def plot(self, **kwargs):
        """Plot xvg file data.

        The first column of the data is always taken as the abscissa
        X. Additional columns are plotted as ordinates Y1, Y2, ...

        In the special case that there is only a single column then this column
        is plotted against the index, i.e. (N, Y).

        :Arguments:
          - *transform*: function *transform(array) --> array* which transforms
            the original array; must return a 2D numpy array of shape [X,
            Y1, Y2, ...] where X, Y1, ... are column vectors.
            By default the transformation is the identity.
          - All other keyword arguments are passed on to :func:`pylab.plot`.
        """
        import pylab
        transform = kwargs.pop('transform', lambda x: x)  # default is identity transformation
        a = numpy.asarray(transform(self.array))
        if len(a.shape) == 1:
            # special case: plot against index; plot would do this automatically but 
            # we'll just produce our own xdata and pretend that this was X all along
            X = numpy.arange(len(a))
            a = numpy.concatenate([[X], [a]])  # does NOT overwrite original a but make a new one
        kwargs['xdata'] = a[0]          # abscissa set separately
        pylab.plot(a[1:].T, **kwargs)   # plot all other columns in parallel
        

def iterable(obj):
    """Returns ``True`` if *obj* can be iterated over and is *not* a  string."""
    if type(obj) is str:
        return False    # avoid iterating over characters of a string

    if hasattr(obj, 'next'):
        return True    # any iterator will do 
    try: 
        len(obj)       # anything else that might work
    except TypeError: 
        return False
    return True

def asiterable(obj):
    """Returns obj so that it can be iterated over; a string is *not* treated as iterable"""
    if not iterable(obj):
        obj = [obj]
    return obj


# In utilities so that it can be safely used in tools, cbook, ...

def unlink_f(path):
    """Unlink path but do not complain if file does not exist."""
    try:
        os.unlink(path)
    except OSError, err:
        if err.errno <> errno.ENOENT:
            raise

def unlink_gmx(*args):
    """Unlink (remove) Gromacs file(s) and all corresponding backups."""
    for path in args:
        unlink_f(path)
    unlink_gmx_backups(*args)

def unlink_gmx_backups(*args):
    """Unlink (rm) all backup files corresponding to the listed files."""
    for path in args:
        dirname, filename = os.path.split(path)
        fbaks = glob.glob(os.path.join(dirname, '#'+filename+'.*#'))
        for bak in fbaks:
            unlink_f(bak)

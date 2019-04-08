# GromacsWrapper: utilities.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
:mod:`gromacs.utilities` -- Helper functions and classes
========================================================

The module defines some convenience functions and classes that are
used in other modules; they do *not* make use of :mod:`gromacs.tools`
or :mod:`gromacs.cbook` and can be safely imported at any time.


Classes
-------

:class:`FileUtils` provides functions related to filename handling. It
can be used as a base or mixin class. The :class:`gromacs.analysis.Simulation`
class is derived from it.

.. autoclass:: FileUtils
   :members:
.. autoclass:: AttributeDict
.. autoclass:: Timedelta

Functions
---------

Some additional convenience functions that deal with files and
directories:

.. function:: openany(directory[,mode='r'])

   Context manager to open a compressed (bzip2, gzip) or plain file
   (uses :func:`anyopen`).

.. autofunction:: anyopen
.. autofunction:: isstream
.. autofunction:: realpath
.. function:: in_dir(directory[,create=True])

   Context manager to execute a code block in a directory.

   * The *directory* is created if it does not exist (unless
     *create* = ``False`` is set)
   * At the end or after an exception code always returns to
     the directory that was the current directory before entering
     the block.

.. autofunction:: find_first
.. autofunction:: withextsep
.. autofunction:: which

Functions that improve list processing and which do *not* treat
strings as lists:

.. autofunction:: iterable
.. autofunction:: asiterable
.. autofunction:: firstof

Functions that help handling Gromacs files:

.. autofunction:: unlink_f
.. autofunction:: unlink_gmx
.. autofunction:: unlink_gmx_backups
.. autofunction:: number_pdbs

Functions that make working with matplotlib_ easier:

.. _matplotlib: http://matplotlib.sourceforge.net/

.. autofunction:: activate_subplot
.. autofunction:: remove_legend


Miscellaneous functions:

.. autofunction:: convert_aa_code
.. autofunction:: autoconvert
.. autofunction:: hasmethod

Data
----

.. autodata:: amino_acid_codes

"""
from __future__ import absolute_import, with_statement, division

__docformat__ = "restructuredtext en"

import six
from six import string_types

import sys
import os
import glob
import fnmatch
import itertools
import re
import warnings
import errno
import subprocess
from contextlib import contextmanager
import bz2
import gzip
import datetime

import numpy

import logging
logger = logging.getLogger('gromacs.utilities')

from .exceptions import AutoCorrectionWarning


class AttributeDict(dict):
    """A dictionary with pythonic access to keys as attributes --- useful for interactive work."""
    def __getattribute__(self, x):
        try:
            return super(AttributeDict,self).__getattribute__(x)
        except AttributeError:
            return self[x]
    def __setattr__(self, name, value):
        try:
            super(AttributeDict, self).__setitem__(name, value)
        except KeyError:
            super(AttributeDict, self).__setattr__(name, value)

    def __getstate__(self):
        return self

    def __setstate__(self, state):
        self.update(state)


def autoconvert(s):
    """Convert input to a numerical type if possible.

    1. A non-string object is returned as it is
    2. Try conversion to int, float, str.
    """
    if type(s) is not str:
        return s
    for converter in int, float, str:   # try them in increasing order of lenience
        try:
            s = [converter(i) for i in s.split()]
            if len(s) == 1:
                return s[0]
            else:
                return numpy.array(s)
        except (ValueError, AttributeError):
            pass
    raise ValueError("Failed to autoconvert {0!r}".format(s))

@contextmanager
def openany(datasource, mode='rt', reset=True):
    """Context manager for :func:`anyopen`.

    Open the `datasource` and close it when the context of the :keyword:`with`
    statement exits.

    `datasource` can be a filename or a stream (see :func:`isstream`). A stream
    is reset to its start if possible (via :meth:`~io.IOBase.seek` or
    :meth:`~cString.StringIO.reset`).

    The advantage of this function is that very different input sources
    ("streams") can be used for a "file", ranging from files on disk (including
    compressed files) to open file objects to sockets and strings---as long as
    they have a file-like interface.

    :Arguments:
      *datasource*
           a file or a stream
      *mode*
           {'r', 'w'} (optional), open in r(ead) or w(rite) mode
      *reset*
           bool (optional) try to read (`mode` 'r') the stream from the
           start [``True``]


    **Example**

    Open a gzipped file and process it line by line::

        with openany("input.pdb.gz") as pdb:
            for line in pdb:
                if line.startswith('ATOM'):
                    print(line)

    Open a URL and read it::

       import urllib2
       with openany(urllib2.urlopen("https://www.mdanalysis.org/")) as html:
           print(html.read())


    .. SeeAlso::
       :func:`anyopen`

    """
    stream = anyopen(datasource, mode=mode, reset=reset)
    try:
        yield stream
    finally:
        stream.close()


# On python 3, we want to use bz2.open to open and uncompress bz2 files. That
# function allows to specify the type of the uncompressed file (bytes ot text).
# The function does not exist in python 2, thus we must use bz2.BZFile to
# which we cannot tell if the uncompressed file contains bytes or text.
# Therefore, on python 2 we use a proxy function that removes the type of the
# uncompressed file from the `mode` argument.
try:
    bz2.open
except AttributeError:
    # We are on python 2 and bz2.open is not available
    def bz2_open(filename, mode):
        """Open and uncompress a BZ2 file"""
        mode = mode.replace('t', '').replace('b', '')
        return bz2.BZ2File(filename, mode)
else:
    # We are on python 3 so we can use bz2.open
    bz2_open = bz2.open


def anyopen(datasource, mode='rt', reset=True):
    """Open datasource (gzipped, bzipped, uncompressed) and return a stream.

    `datasource` can be a filename or a stream (see :func:`isstream`). By
    default, a stream is reset to its start if possible (via
    :meth:`~io.IOBase.seek` or :meth:`~cString.StringIO.reset`).

    If possible, the attribute ``stream.name`` is set to the filename or
    "<stream>" if no filename could be associated with the *datasource*.

    :Arguments:
      *datasource*
        a file (from :class:`file` or :func:`open`) or a stream (e.g. from
        :func:`urllib2.urlopen` or :class:`cStringIO.StringIO`)
      *mode*
        {'r', 'w', 'a'} (optional),
        Open in r(ead), w(rite) or a(ppen) mode. More complicated
        modes ('r+', 'w+', ...) are not supported; only the first letter of
        `mode` is used and thus any additional modifiers are silently ignored.
      *reset*
        bool (optional),
        try to read (`mode` 'r') the stream from the start

    :Returns:
       file-like object


    .. SeeAlso::
       :func:`openany` to be used with the :keyword:`with` statement.

    """
    handlers = {'bz2': bz2_open, 'gz': gzip.open, '': open}

    if mode.startswith('r'):
        if isstream(datasource):
            stream = datasource
            try:
                filename = str(stream.name)  # maybe that does not always work?
            except AttributeError:
                filename = "<stream>"
            if reset:
                try:
                    stream.reset()
                except (AttributeError, IOError):
                    try:
                        stream.seek(0)
                    except (AttributeError, IOError):
                        warnings.warn("Stream {0}: not guaranteed to be at the beginning."
                                      "".format(filename),
                                      category=StreamWarning)
        else:
            stream = None
            filename = datasource
            for ext in ('bz2', 'gz', ''):  # file == '' should be last
                openfunc = handlers[ext]
                stream = _get_stream(datasource, openfunc, mode=mode)
                if stream is not None:
                    break
            if stream is None:
                raise IOError(errno.EIO, "Cannot open file or stream in mode={mode!r}.".format(**vars()), repr(filename))
    elif mode.startswith('w') or mode.startswith('a'):  # append 'a' not tested...
        if isstream(datasource):
            stream = datasource
            try:
                filename = str(stream.name)  # maybe that does not always work?
            except AttributeError:
                filename = "<stream>"
        else:
            stream = None
            filename = datasource
            name, ext = os.path.splitext(filename)
            if ext.startswith('.'):
                ext = ext[1:]
            if not ext in ('bz2', 'gz'):
                ext = ''  # anything else but bz2 or gz is just a normal file
            openfunc = handlers[ext]
            stream = openfunc(datasource, mode=mode)
            if stream is None:
                raise IOError(errno.EIO, "Cannot open file or stream in mode={mode!r}.".format(**vars()), repr(filename))
    else:
        raise NotImplementedError("Sorry, mode={mode!r} is not implemented for {datasource!r}".format(**vars()))
    try:
        stream.name = filename
    except (AttributeError, TypeError):
        pass  # can't set name (e.g. cStringIO.StringIO)
    return stream


def _get_stream(filename, openfunction=open, mode='r'):
    """Return open stream if *filename* can be opened with *openfunction* or else ``None``."""
    try:
        stream = openfunction(filename, mode=mode)
    except (IOError, OSError) as err:
        # An exception might be raised due to two reasons, first the openfunction is unable to open the file, in this
        # case we have to ignore the error and return None. Second is when openfunction can't open the file because
        # either the file isn't there or the permissions don't allow access.
        if errno.errorcode[err.errno] in ['ENOENT', 'EACCES']:
            six.reraise(*sys.exc_info())
        return None
    if mode.startswith('r'):
        # additional check for reading (eg can we uncompress) --- is this needed?
        try:
            stream.readline()
        except IOError:
            stream.close()
            stream = None
        except:
            stream.close()
            raise
        else:
            stream.close()
            stream = openfunction(filename, mode=mode)
    return stream

def hasmethod(obj, m):
    """Return ``True`` if object *obj* contains the method *m*.

    .. versionadded:: 0.7.1
    """
    return hasattr(obj, m) and callable(getattr(obj, m))

def isstream(obj):
    """Detect if `obj` is a stream.

    We consider anything a stream that has the methods

    - ``close()``

    and either set of the following

    - ``read()``, ``readline()``, ``readlines()``
    - ``write()``, ``writeline()``, ``writelines()``

    :Arguments:
      *obj*
          stream or str

    :Returns:
      *bool*, ``True`` if `obj` is a stream, ``False`` otherwise

    .. SeeAlso::
       :mod:`io`


    .. versionadded:: 0.7.1
    """
    signature_methods = ("close",)
    alternative_methods = (
        ("read", "readline", "readlines"),
        ("write", "writeline", "writelines"))

    # Must have ALL the signature methods
    for m in signature_methods:
        if not hasmethod(obj, m):
            return False
    # Must have at least one complete set of alternative_methods
    alternative_results = [
        numpy.all([hasmethod(obj, m) for m in alternatives])
        for alternatives in alternative_methods]
    return numpy.any(alternative_results)

# TODO: make it work for non-default charge state amino acids.
#: translation table for 1-letter codes --> 3-letter codes
#: .. Note: This does not work for HISB and non-default charge state aa!
amino_acid_codes = {'A':'ALA', 'C':'CYS', 'D':'ASP', 'E':'GLU',
                    'F':'PHE', 'G':'GLY', 'H':'HIS', 'I':'ILE',
                    'K':'LYS', 'L':'LEU', 'M':'MET', 'N':'ASN',
                    'P':'PRO', 'Q':'GLN', 'R':'ARG', 'S':'SER',
                    'T':'THR', 'V':'VAL', 'W':'TRP', 'Y':'TYR'}
inverse_aa_codes = {three: one for one,three in amino_acid_codes.items()}

def convert_aa_code(x):
    """Converts between 3-letter and 1-letter amino acid codes."""
    if len(x) == 1:
        return amino_acid_codes[x.upper()]
    elif len(x) == 3:
        return inverse_aa_codes[x.upper()]
    else:
        raise ValueError("Can only convert 1-letter or 3-letter amino acid codes, "
                         "not %r" % x)

@contextmanager
def in_dir(directory, create=True):
    """Context manager to execute a code block in a directory.

    * The directory is created if it does not exist (unless
      create=False is set)
    * At the end or after an exception code always returns to
      the directory that was the current directory before entering
      the block.
    """
    startdir = os.getcwd()
    try:
        try:
            os.chdir(directory)
            logger.debug("Working in {directory!r}...".format(**vars()))
        except OSError as err:
            if create and err.errno == errno.ENOENT:
                os.makedirs(directory)
                os.chdir(directory)
                logger.info("Working in {directory!r} (newly created)...".format(**vars()))
            else:
                logger.exception("Failed to start working in {directory!r}.".format(**vars()))
                raise
        yield os.getcwd()
    finally:
        os.chdir(startdir)

def realpath(*args):
    """Join all args and return the real path, rooted at /.

    Expands ``~`` and environment variables such as :envvar:`$HOME`.

    Returns ``None`` if any of the args is none.
    """
    if None in args:
        return None
    return os.path.realpath(
        os.path.expandvars(os.path.expanduser(os.path.join(*args))))

def find_first(filename, suffices=None):
    """Find first *filename* with a suffix from *suffices*.

    :Arguments:
      *filename*
         base filename; this file name is checked first
      *suffices*
         list of suffices that are tried in turn on the root of *filename*; can contain the
         ext separator (:data:`os.path.extsep`) or not

    :Returns: The first match or ``None``.
    """
    # struct is not reliable as it depends on qscript so now we just try everything...

    root,extension = os.path.splitext(filename)
    if suffices is None:
        suffices = []
    else:
        suffices = withextsep(suffices)
    extensions = [extension] + suffices  # native name is first
    for ext in extensions:
        fn = root + ext
        if os.path.exists(fn):
            return fn
    return None

def withextsep(extensions):
    """Return list in which each element is guaranteed to start with :data:`os.path.extsep`."""
    def dottify(x):
        if x.startswith(os.path.extsep):
            return x
        return os.path.extsep + x
    return [dottify(x) for x in asiterable(extensions)]

def find_files(directory, pattern):
    """Find files recursively under *directory*, matching *pattern* (generator).

    *pattern* is a UNIX-style glob pattern as used ny :func:`fnmatch.fnmatch`.

    Recipe by Bruno Oliveira from
    http://stackoverflow.com/questions/2186525/use-a-glob-to-find-files-recursively-in-python
    """
    for root, dirs, files in os.walk(directory):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                yield filename

def which(program):
    """Determine full path of executable *program* on :envvar:`PATH`.

    (Jay at http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python)

    .. versionadded:: 0.5.1
    """

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        real_program = realpath(program)
        if is_exe(real_program):
            return real_program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

class FileUtils(object):
    """Mixin class to provide additional file-related capabilities."""

    #: Default extension for files read/written by this class.
    default_extension = None

    def _init_filename(self, filename=None, ext=None):
        """Initialize the current filename :attr:`FileUtils.real_filename` of the object.

        Bit of a hack.

        - The first invocation must have ``filename != None``; this will set a
          default filename with suffix :attr:`FileUtils.default_extension`
          unless another one was supplied.

        - Subsequent invocations either change the filename accordingly or
          ensure that the default filename is set with the proper suffix.

        """

        extension = ext or self.default_extension
        filename = self.filename(filename, ext=extension, use_my_ext=True, set_default=True)
        #: Current full path of the object for reading and writing I/O.
        self.real_filename = os.path.realpath(filename)

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

        .. versionchanged:: 0.3.1
           An empty string as *ext* = "" will suppress appending an extension.
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
            if ext.startswith(os.extsep):
                ext = ext[1:]  # strip a dot to avoid annoying mistakes
            if ext != "":
                filename = filename + os.extsep + ext
        return filename

    def check_file_exists(self, filename, resolve='exception', force=None):
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

        Alternatively, set *force* for the following behaviour (which
        ignores *resolve*):

        ``True``
              same as *resolve* = "ignore" (will allow overwriting of files)
        ``False``
              same as *resolve* = "exception" (will prevent overwriting of files)
        ``None``
              ignored, do whatever *resolve* says
        """
        def _warn(x):
            msg = "File {0!r} already exists.".format(x)
            logger.warn(msg)
            warnings.warn(msg)
            return True
        def _raise(x):
            msg = "File {0!r} already exists.".format(x)
            logger.error(msg)
            raise IOError(errno.EEXIST, x, msg)
        solutions = {'ignore': lambda x: False,      # file exists, but we pretend that it doesn't
                     'indicate': lambda x: True,     # yes, file exists
                     'warn': _warn,
                     'warning': _warn,
                     'exception': _raise,
                     'raise': _raise,
                     }

        if force is True:
            resolve = 'ignore'
        elif force is False:
            resolve = 'exception'

        if not os.path.isfile(filename):
            return False
        else:
            return solutions[resolve](filename)

    def infix_filename(self, name, default, infix, ext=None):
        """Unless *name* is provided, insert *infix* before the extension *ext* of *default*."""
        if name is None:
            p, oldext = os.path.splitext(default)
            if ext is None:
                ext = oldext
            if ext.startswith(os.extsep):
                ext = ext[1:]
            name = self.filename(p+infix, ext=ext)
        return name

    def __repr__(self):
        fmt = "{0!s}(filename=%r)".format(self.__class__.__name__)
        try:
            fn =  self.filename()
        except ValueError:
            fn = None
        return fmt % fn


def iterable(obj):
    """Returns ``True`` if *obj* can be iterated over and is *not* a  string."""
    if isinstance(obj, string_types):
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

def firstof(obj):
    """Returns the first entry of a sequence or the obj.

    Treats strings as single objects.
    """
    return asiterable(obj)[0]

# In utilities so that it can be safely used in tools, cbook, ...

def unlink_f(path):
    """Unlink path but do not complain if file does not exist."""
    try:
        os.unlink(path)
    except OSError as err:
        if err.errno != errno.ENOENT:
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

def mkdir_p(path):
    """Create a directory *path* with subdirs but do not complain if it exists.

    This is like GNU ``mkdir -p path``.
    """
    try:
        os.makedirs(path)
    except OSError as err:
        if err.errno != errno.EEXIST:
            raise

def cat(f=None, o=None):
    """Concatenate files *f*=[...] and write to *o*"""
    # need f, o to be compatible with trjcat and eneconv
    if f is None or o is None:
        return
    target = o
    infiles = asiterable(f)
    logger.debug("cat {0!s} > {1!s} ".format(" ".join(infiles), target))
    with open(target, 'w') as out:
        rc = subprocess.call(['cat'] + infiles, stdout=out)
    if rc != 0:
        msg = "failed with return code {0:d}: cat {1!r} > {2!r} ".format(rc, " ".join(infiles), target)
        logger.exception(msg)
        raise OSError(errno.EIO, msg, target)


# helpers for matplotlib
def activate_subplot(numPlot):
    """Make subplot *numPlot* active on the canvas.

    Use this if a simple ``subplot(numRows, numCols, numPlot)``
    overwrites the subplot instead of activating it.
    """
    # see http://www.mail-archive.com/matplotlib-users@lists.sourceforge.net/msg07156.html
    from pylab import gcf, axes
    numPlot -= 1  # index is 0-based, plots are 1-based
    return axes(gcf().get_axes()[numPlot])

def remove_legend(ax=None):
    """Remove legend for axes or gca.

    See http://osdir.com/ml/python.matplotlib.general/2005-07/msg00285.html
    """
    from pylab import gca, draw
    if ax is None:
        ax = gca()
    ax.legend_ = None
    draw()


# time functions
class Timedelta(datetime.timedelta):
    """Extension of :class:`datetime.timedelta`.

    Provides attributes ddays, dhours, dminutes, dseconds to measure
    the delta in normal time units.

    ashours gives the total time in fractional hours.
    """

    @property
    def dhours(self):
        """Hours component of the timedelta."""
        return self.seconds // 3600

    @property
    def dminutes(self):
        """Minutes component of the timedelta."""
        return self.seconds // 60 - 60*self.dhours

    @property
    def dseconds(self):
        """Seconds component of the timedelta."""
        return self.seconds - 3600*self.dhours - 60*self.dminutes

    @property
    def ashours(self):
        """Timedelta in (fractional) hours."""
        return 24*self.days + self.seconds / 3600.0

    def strftime(self, fmt="%d:%H:%M:%S"):
        """Primitive string formatter.

        The only directives understood are the following:
          ============   ==========================
          Directive      meaning
          ============   ==========================
          %d             day as integer
          %H             hour  [00-23]
          %h             hours including days
          %M             minute as integer [00-59]
          %S             second as integer [00-59]
          ============   ==========================
        """
        substitutions = {
            "%d": str(self.days),
            "%H": "{0:02d}".format(self.dhours),
            "%h": str(24*self.days + self.dhours),
            "%M": "{0:02d}".format(self.dminutes),
            "%S": "{0:02d}".format(self.dseconds),
            }
        s = fmt
        for search, replacement in substitutions.items():
            s = s.replace(search, replacement)
        return s


NUMBERED_PDB = re.compile(r"(?P<PREFIX>.*\D)(?P<NUMBER>\d+)\.(?P<SUFFIX>pdb)")

def number_pdbs(*args, **kwargs):
    """Rename pdbs x1.pdb ... x345.pdb --> x0001.pdb ... x0345.pdb

    :Arguments:
       - *args*: filenames or glob patterns (such as "pdb/md*.pdb")
       - *format*: format string including keyword *num* ["%(num)04d"]
    """

    format = kwargs.pop('format', "%(num)04d")
    name_format = "%(prefix)s" + format +".%(suffix)s"

    for f in itertools.chain.from_iterable(map(glob.glob, args)):
        m = NUMBERED_PDB.search(f)
        if m is None:
            continue
        num = int(m.group('NUMBER'))
        prefix = m.group('PREFIX')
        suffix = m.group('SUFFIX')
        newname = name_format % vars()
        logger.info("Renaming {f!r} --> {newname!r}".format(**vars()))
        try:
            os.rename(f, newname)
        except OSError:
            logger.exception("renaming failed")

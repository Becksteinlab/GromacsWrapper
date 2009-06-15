# $Id$
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""Helper functions."""

import os.path
import warnings
import bz2, gzip

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


class FileUtils(object):
    """Mixin class to provide additional IO capabilities."""

    def filename(self,filename=None,ext=None,set_default=False,use_my_ext=False):
        """Supply a file name for the object.

        fn = filename()             ---> <default_filename>
        fn = filename('name.ext')   ---> 'name'
        fn = filename(ext='pickle') ---> <default_filename>'.pickle'
        fn = filename('name.inp','pdf') --> 'name.pdf'
        fn = filename('foo.pdf',ext='png',use_my_ext=True) --> 'foo.pdf'

        The returned filename is stripped of the extension (use_my_ext=False) and
        if provided, another extension is appended. Chooses a default if no
        filename is given.  Raises a ValueError exception if no default file name
        is known.

        If set_default=True then the default filename is also set.

        use_my_ext=True lets the suffix of a provided filename take priority over a
        default ext(tension).
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
        """If a file exists then continue with the recipe specified in resolve.

        resolve must be one of
           ignore:           always return False
           indicate:         return True if it exists
           warn:             indicate and issue a UserWarning
           exception:        raise IOError if it exists
        """
        def _warn(x):
            warnings.warn("File %r already exists." % x)
            return True
        def _raise(x):
            raise IOError("File %r already exists." % x)
        solutions = {'ignore': lambda x: False,      # file exists, but we pretend that it doesn't
                     'indicate': lambda x: True,     # yes, file exists
                     'warn': _warn,
                     'exception': _raise,
                     }
        if not os.path.isfile(filename):
            return False
        else:
            return solutions[resolve](filename)


def iterable(obj):
    """Returns True if obj can be iterated over and is NOT a  string."""
    if type(obj) is str:
        return False    # avoid iterating over characters of a string

    if hasattr(obj, 'next'):
        return True    # any iterator will do 
    try: 
        len(obj)       # anything else that might work
    except TypeError: 
        return False
    return True

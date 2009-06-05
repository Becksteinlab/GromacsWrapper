# $Id$
"""Helper functions."""

import bz2, gzip

def anyopen(datasource):
    """Open datasource and return a stream."""
    if hasattr(datasource,'next') or hasattr(datasource,'readline'):
        stream = datasource
        filename = '(%s)' % stream.name  # maybe that does not always work?        
    else:
        stream = None
        filename = datasource
        for openfunc in bz2.BZ2File, gzip.open, file:
            stream = _trialreadline(datasource)
            if not stream is None:
                break
        if stream is None:
            raise IOError("Cannot open %r for reading." % filename)
    return stream, filename

def _trialreadline(filename, openfunction=file):
    try:
        stream = openfunction(filename,'r')
        stream.readline()
        stream.close()
        stream = bz2.BZ2File(self.filename,'r')
    except IOError:
        stream.close()
        stream = None
    return stream

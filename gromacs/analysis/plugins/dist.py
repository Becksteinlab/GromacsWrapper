# $Id$
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
``analysis.plugins.dist`` --- Helper Class for ``g_dist``
=========================================================

:mod:`dist` contains helper classes for other analysis plugins that
want to make use of the Gromacs command ``g_dist``.


Overview
--------

The task we are solving is to analyze output from ::

   g_dist -f md.xtc -s md.tpr -n cys_ow.ndx -dist 1.0 | bzip2 -vc > mindist_C60_OW_1nm.dat.bz2

and produce a histogram of minimum contact distances. This should
provide an estimate for water accessibility of the atom (here: SG of
Cys60).

File format
-----------

``g_dist`` with the ``-dist CUTOFF`` option writes to stdout the
identity of all atoms within the cutoff distance and the distance
itself::

   Selected 22: 'CYSH_CYSH_60_&_SG'
   Selected 25: 'OW'
   ....
   t: 184  6682 SOL 35993 OW  0.955138 (nm)
   t: 184  10028 SOL 46031 OW  0.803889 (nm)
   t: 185  6682 SOL 35993 OW  0.879949 (nm)
   t: 185  10028 SOL 46031 OW  0.738299 (nm)
   t: 186  6682 SOL 35993 OW  0.897016 (nm)
   t: 186  10028 SOL 46031 OW  0.788268 (nm)
   t: 187  6682 SOL 35993 OW  0.997688 (nm)
   ....


Classes
-------

.. autoclass:: Mindist
   :members: histogram, hist, dist, edges, midpoints, plot
.. autoclass:: GdistData
   :members: __iter__

"""
__docformat__ = "restructuredtext en"

import re
import numpy

from recsql import SQLarray   # my own ReqSQL module

import gromacs.utilities

class Mindist(object):
    """The Mindist class allows analysis of the output from ``g_dist -dist CUTOFF``.

    Output is read from a file or stream. The raw data is transformed
    into a true 'mindist' time series (available in the
    :attr:`Mindist.distances` attribute): for each frame only the
    shortest distance is stored (whereas ``g_dist`` provides *all*
    distances below the cutoff).

    :TODO:
      * Save analysis to pickle or data files.
      * Export data as simple data files for plotting in other programs.

    .. Note::

       :class:`gromacs.tools.G_mindist` is apparently providing
       exactly the service that is required: a timeseries of the
       minimum distance between two groups. Feel free to use that tool
       instead :-).

    """

    def __init__(self,datasource,cutoff=None):
        """Read mindist data from file or stream.

        :Arguments:
          *datasource*
             a filename (plain, gzip, bzip2) or file object
          *cutoff*
             the ``-dist CUTOFF`` that was provided to ``g_dist``; if supplied
             we work around a bug in ``g_dist`` (at least in Gromacs 4.0.2) in which
             sometimes numbers >> CUTOFF are printed.
        """
        stream, self.filename = gromacs.utilities.anyopen(datasource)
        try:
            M = GdistData(stream)
            # BIG database in memory ... can be accessed via SQL
            all_distances = SQLarray('distances', records=M, columns=('frame','distance'))
        finally:
            stream.close()
        if cutoff is None:
            cutoff_filter = ""
        else:
            cutoff_filter = "WHERE distance <= %d" % float(cutoff)
        self.distances = all_distances.selection(
            "SELECT frame, MIN(distance) AS distance FROM __self__ "+cutoff_filter+" GROUP BY frame",
            name="mindistances", cache=False)

    def histogram(self,nbins=None,lo=None,hi=None,midpoints=False,normed=True):
        """Returns a distribution or histogram of the minimum distances.

        If no values for the bin edges are given then they are set to
        0.1 below and 0.1 above the minimum and maximum values seen in
        the data.

        If the number of bins is not provided then it is set so that
        on average 100 counts come to a bin. Set nbins manually if the
        histogram only contains a single bin (and then get more data)!

        :Keywords:

           *nbins* : int
              number of bins
           *lo* : float
              lower edge of histogram
           *hi* : float
              upper edge of histogram
           *midpoints* : boolean
              False: return edges. True: return midpoints
           *normed* : boolean
              True: return probability distribution. False: histogram

        """
        D = self.distances
        if lo is None or hi is None:
            dmin,dmax = D.limits('distance')
            if lo is None:
                lo = round(dmin - 0.1, 1)
                if lo < 0:
                    lo = 0.0
            if hi is None:
                hi = round(dmax + 0.1, 1)
        if nbins is None:
            nbins = int(len(D)/100)
        FUNC = 'distribution'
        if not normed:
            FUNC = 'histogram'
        SQL = """SELECT %(FUNC)s(distance,%(nbins)d,%(lo)f,%(hi)f) AS "h [Object]"
                 FROM __self__""" % vars()
        (((h,e),),) = D.sql(SQL, asrecarray=False)
        # should probably cache this...
        if normed:
            self.__distribution = h
        else:
            self.__histogram = h
        self.__edges = e
        if midpoints:
            return h, self.midpoints
        return h, e

    def hist():
        doc = "Histogram of the minimum distances."
        def fget(self):
            try:
                return self.__histogram
            except AttributeError:
                raise RuntimeError("run histogram(...,normed=False) first")
        return locals()
    hist = property(**hist())

    def dist():
        doc = "Distribution of the minimum distances."
        def fget(self):
            try:
                return self.__distribution
            except AttributeError:
                raise RuntimeError("run histogram(...,normed=True) first")
        return locals()
    dist = property(**dist())

    def edges():
        doc = "Edges of the histogram of the minimum distances."
        def fget(self):
            try:
                return self.__edges
            except AttributeError:
                raise RuntimeError("run histogram(...) first")
        return locals()
    edges = property(**edges())

    def midpoints():
        doc = "Midpoints of the histogram of the minimum distances."
        def fget(self):
            e = self.__edges
            return 0.5*(e[:-1] + e[1:])
        return locals()
    midpoints = property(**midpoints())

    def plot(self,**kwargs):
        """Plot histograms with matplotlib's plot() function::

           plot(**histogramargs, **plotargs)

        Arguments for both :meth:`Mindist.histogram` and :func:`pylab.plot` can be provided (qv).
        """
        import pylab
        kwargs['midpoints'] = True
        histogram_args = ('nbins','lo','hi','midpoints','normed')
        histargs = kwargs.copy()
        histargs = {k: kwargs.pop(k) for k in histargs if k in histogram_args}
        h,m = self.histogram(**histargs)

        pylab.plot(m, h, **kwargs)
        pylab.xlabel('minimum distance (nm)')

class GdistData(object):
    """Object that represents the standard output of ``g_dist -dist CUTOFF``.

    Initialize from a stream (e.g. a pipe) and the iterate over the
    instance to get the data line by line. Each line consists of a
    tuple

      (*frame*, *distance*)
    """
    data_pattern = re.compile("""t:\s*                 # marker for beginning of line
                (?P<FRAME>\d+)\s+                      # frame number (?)
                (?P<RESID>\d+)\s+(?P<RESNAME>\w+)\s+   # resid and residue name
                (?P<ATOMID>\d+)\s+(?P<ATOMNAME>\w+)\s+ # atomid and atom name
                (?P<DISTANCE>[0-9]+\.?[0-9]+)          # distance
                \s+\(nm\)""", re.VERBOSE)

    def __init__(self,stream):
        """Initialize with an open stream to the data (eg stdin or file).

        :Arguments:
           *stream*
              open stream (file or pipe or really any iterator
              providing the data from ``g_dist``); the *stream* is not
              closed automatically when the iterator completes.
        """
        self.stream = stream

    def __iter__(self):
        """Iterator that filters the input stream and returns (frame, distance) tuples."""
        for line in self.stream:
            line = line.strip()
            m = self.data_pattern.search(line)
            if m:
                distance = float(m.group('DISTANCE'))
                frame = int(m.group('FRAME'))
                yield frame,distance

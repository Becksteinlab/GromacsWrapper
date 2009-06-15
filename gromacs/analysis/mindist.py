# $Id$
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
Overview
========

Analyze output from::

 printf '22\n25\n' | \
   g_dist -f ../md.xtc -s ../md.tpr -n cys_ow.ndx -dist 1.0 | \
   bzip2 -vc > mindist_C60_OW_1nm.dat.bz2 

and produce a histogram of minimum contact distances. This should
provide an estimate for water accessibility of the atom (here: SG of
Cys60).

File format
===========
  
``g_mindist`` with the ``-dist CUTOFF`` option writes to stdout the
identity of all atoms within the cutoff distance and the distance
itself.::

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
   ...

"""

import re
import numpy

from recsql import SQLarray   # my own ReqSQL module
from recsql.sqlutil import FakeRecArray

import gromacs.utilities

class Mindist(object):
    """The Mindist class allows analysis of the output from ``g_dist -dist CUTOFF``.

    Output is read from a file or stream. The raw data (attribute
    ``all_distances``) is transformed into a true 'mindist' time
    series (available in the ``distances`` attribute): for each frame
    only the shortest distance is stored (whereas g_mindist provides
    *all* distances below the cutoff).

    :Attributes:
    all_distances      data from g_mindist (frame, distance) 
    distances          time (frame) series of the shortest distances

    :Methods:
    histogram          histogram of the mindist time series
    plot               compute histograms and plot with matplotlib

    :TODO:
    * Save analysis to pickle or data files.
    * Export data as simple data files for plotting in other programs.
    * Make implementation more memory efficient.
    """

    def __init__(self,datasource,cutoff=None):
        """Read mindist data from file or stream.

        M = Mindist(datasource, cutoff=1.0)
        
        :Arguments:
        datasource:     a filename (plain, gzip, bzip2) or file object
        cutoff:         the ``-dist CUTOFF`` that was provided to ``g_dist``; if supplied 
                        we work around a bug in ``g_dist`` (at least in Gromacs 4.0.2) in which
                         sometimes numbers >> CUTOFF are printed.
        """
        stream, self.filename = gromacs.utilities.anyopen(datasource)
        try:
            M = GdistData(stream)
            # make a tmp recarray (BIG array in memory!!) to pass into
            # db (crappy...)
            #_distances = numpy.rec.fromrecords(
            #    [(frame,distance) for frame,distance in M],  # pull data from file
            #    names='frame,distance', formats='i4,f8')     # MUST be i4 (i8 gives sqlite unsuported type error)
            _distances = FakeRecArray(M, columns=('frame','distance'))
            # BIG database in memory!!
            all_distances = SQLarray('distances', _distances)  # ... can be accessed via SQL
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

        hist,edges = histogram(nbins=10, hi=1)

        If no values for the bin edges are given then they are set to
        0.1 below and 0.1 above the minimum and maximum values seen in
        the data.

        If the number of bins is not provided then it is set so that
        on average 100 counts come to a bin. Set nbins manually if the
        histogram only contains a single bin (and then get more data)!

        :Keyword arguments:
        nbins       number of bins
        lo          lower edge of histogram
        hi          upper edge of histogram
        midpoints   False: return edges. True: return midpoints
        normed      True: return probability distribution. False: histogram
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
        """Plot histograms with matplotlib's plot() function.

        plot(**histogramargs, **plotargs)

        Arguments for both histogram() and plot() can be provided (qv).
        """
        import pylab
        kwargs['midpoints'] = True
        histogram_args = ('nbins','lo','hi','midpoints','normed')
        histargs = kwargs.copy()
        histargs = dict([(k, kwargs.pop(k)) for k in histargs if k in histogram_args])
        h,m = self.histogram(**histargs)        

        pylab.plot(m, h, **kwargs)
        pylab.xlabel('minimum distance (nm)')

class GdistData(object):
    """Object that represents the output of g_dist -dist CUTOFF"""
    data_pattern = re.compile("""t:\s*                 # marker for beginning of line
                (?P<FRAME>\d+)\s+                      # frame number (?)
                (?P<RESID>\d+)\s+(?P<RESNAME>\w+)\s+   # resid and residue name
                (?P<ATOMID>\d+)\s+(?P<ATOMNAME>\w+)\s+ # atomid and atom name
                (?P<DISTANCE>[0-9]+\.?[0-9]+)          # distance
                \s+\(nm\)""", re.VERBOSE)

    def __init__(self,stream):
        """Initialize with an open stream to the data (eg stdin or file)"""
        self.stream = stream

    def __iter__(self):
        for line in self.stream:
            line = line.strip()
            m = self.data_pattern.search(line)
            if m:
                distance = float(m.group('DISTANCE'))
                frame = int(m.group('FRAME'))
                yield frame,distance

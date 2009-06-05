#!/usr/bin/env python
# $Id$
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

import bz2
import re
import numpy

from recsql import SQLarray   # my own ReqSQL module

import gromacs.utilities

class Mindist(object):
    """The Mindist class allows analysis of the output from ``g_dist -dist CUTOFF``.

    Output is read from a bzipped file. The raw data (attribute
    ``all_distances``) is transformed into a true 'mindist' time
    series (available in the ``distances`` attribute): for each frame
    only the shortest distance is stored (whereas g_mindist provides
    *all* distances below the cutoff).

    :Attributes:
    all_distances      data from g_mindist (frame, distance) 
    distances          time (frame) series of the shortest distances

    :Methods:
    histogram          histogram of the mindist time series
    """

    def __init__(self,datasource):
        """Read mindist data from file or stream."""

        self.filename, stream = gromacs.utilities.anyopen(datasource)
        M = GdistData(stream)
        _distances = numpy.rec.fromrecords(
            [(frame,distance) for frame,distance in M],  # pull data from file
            names='frame,distance')                      # and make a tmp recarray
        self.all_distances = SQLarray('distances', _distances)  # ... can be accessed via SQL
        self.distances = self.all_distances.selection(
            'SELECT frame, MIN(distance) AS distance FROM __self__ GROUP BY frame',
            name="mindistances")
        
    def histogram(self,nbins=None,lo=None,hi=None,midpoints=False,normed=True):
        """Returns a histogram of the minimum distances.

        hist,edges = histogram(nbins=10, hi=1)

        If no values for the bin edges are given then they are set to
        0.1 below and 0.1 above the minimum and maximum values seen in
        the data.

        If the number of bins is not provided then it is set so that
        on average 100 counts come to a bin.

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
        if midpoints:
            e = 0.5*(e[:-1] + e[1:])
        return h,e

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

# $Id$
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
:mod:`gromacs.analysis.plugins.gridmatmd` --- Lipid bilayer analysis helper
===========================================================================

This helper module contains code to drive ``GridMAT-MD.pl``, available
from the `GridMAT-MD`_ home page and written by WJ Allen et al
[Allen2009]_ . The GromacsWrapper distribution comes with version
1.0.2 of ``GridMAT-MD.pl`` and includes a small patch so that it can
accept filenames on the command line.



References
----------

.. [Allen2009]  W. J. Allen, J. A. Lemkul, and D. R. Bevan. (2009) "GridMAT-MD: A
                Grid-based Membrane Analysis Tool for Use With Molecular Dynamics."
                J. Comput. Chem. 30 (12): 1952-1958.

.. _GridMAT-MD: http://www.bevanlab.biochem.vt.edu/GridMAT-MD/index.html

Module contents
---------------

.. autoclass:: GridMatMD
   :members:

.. autoclass:: GridMatData
   :members:

"""

import sys
import subprocess
import re
import numpy

import gromacs.tools

class GridMatMD(object):    
    """Analysis of lipid bilayers with GridMAT-MD.

    It requires a configuration file and a list of structure files
    (gro or pdb) as input. See the documentation (pdf_) for the format of the
    config file. Note that the *bilayer* keyword will be ignored in
    the config file.

    .. pdf: http://www.bevanlab.biochem.vt.edu/GridMAT-MD/doc/GridMAT-MD_ug_v1.0.2.pdf
    """
    
    # regular expression to pull useful information from output
    pNX = re.compile("There are (?P<NX>\d+) grid points in the X direction, "
                     "spaced every (?P<DX>[0-9.]+) nanometers")
    pNY = re.compile("There are (?P<NY>\d+) grid points in the Y direction, "
                     "spaced every (?P<DY>[0-9.]+) nanometers")
    pTOP = re.compile('The top leaflet "thickness" will be printed to (?P<FILENAME>.*_top_.*\.dat)')
    pAVG = re.compile('The average bilayer "thickness" will be printed to (?P<FILENAME>.*_average_.*\.dat)')
    pBOT = re.compile('The bottom leaflet "thickness" will be printed to (?P<FILENAME>.*_bottom_.*\.dat)')


    def __init__(self, config, filenames):
        self.config = config
        self.filenames = filenames   # list of filenames (gro or pdb)
        #: Number of bins in the *x* direction. 
        self.nx = None
        #: Numberof bins in the *y* direction.         
        self.ny = None
        #: Bin width in *x*.
        self.dx = None
        #: Bin width in *y*.
        self.dy = None
        #: thickness results (dict of :class:`GridMatData` instances)
        self.thickness = {}
        #: averages of thickness, indexed by "top", "bottom", "average"
        self.arrays = {}
        #: averages of bins
        self.bins = {}
        #: averages of midpoints
        self.midpoints = {}

        self.GridMatMD = gromacs.tools.GridMAT_MD(self.config)

    def run_frame(self, frame):
        """Run GridMAT-MD on a single *frame* and store results."""
        
        thickness = {}
        print "Analyzing frame %r (takes a few seconds)..." % frame
        sys.stdout.flush()
        rc,output,stderr = self.GridMatMD(frame, stdout=False)
        for line in output.split('\n'):
            print '>>> '+line
            m = self.pNX.match(line)
            if m:
                self.nx, self.dx = int(m.group('NX')), float(m.group('DX'))
                continue
            m = self.pNY.match(line)
            if m:
                self.ny, self.dy = int(m.group('NY')), float(m.group('DY'))
                continue
            m = self.pTOP.match(line)
            if m:
                thickness['top'] = m.group('FILENAME')
                continue
            m = self.pBOT.match(line)
            if m:
                thickness['bottom'] = m.group('FILENAME')
                continue
            m = self.pAVG.match(line)
            if m:
                thickness['average'] = m.group('FILENAME')
                continue

        d = {}  # thickness results
        for name, filename in thickness.items():
            d[name] = GridMatData(filename, shape=(self.nx, self.ny),
                                  delta=(self.dx, self.dy))
        self.thickness = d
        return d

    def run(self):
        """Run analysis on all files and average results."""

        # ugly... accumulate averages and divide at the end
        arrays = {}
        bins = {}
        midpoints = {}
        for f in self.filenames:
            d = self.run_frame(f)        # <-- run GridMAT-MD
            for name, data in d.items():
                try:
                    arrays[name] += data.array
                    for dim in xrange(len(bins[name])):
                        bins[name][dim] += data.bins[dim]
                        midpoints[name][dim] += data.midpoints[dim]
                except KeyError:
                    arrays[name] = data.array
                    bins[name] = data.bins
                    midpoints[name] = data.midpoints
        N = float(len(self.filenames))
        for name in arrays:
            arrays[name] /= N
            for dim in xrange(len(bins[name])):
                bins[name][dim] /= N
                midpoints[name][dim] /= N
        self.arrays = arrays
        self.bins = bins
        self.midpoints = midpoints


class GridMatData(object):
    """Represent GridMatMD data file.

    The loaded array data is accessible as a numpy array in
    :attr:`GridMatData.array` and bins and midpoints as
    :attr:`GridMatData.bins` and :attr:`GridMatData.midpoints`
    respectively.
    """

    DATANAME = re.compile("""
                         (?P<NX>\d+)     # number of grid points in X
                         x
                         (?P<NY>\d+)     # number of grid points in Y
                         _(top|bottom|average)_.*\.dat  # all the rest
                         """, re.VERBOSE)


    def __init__(self, filename, shape=None, delta=None):
        """Load the data into a numpy array.

        The *filename* is an output file from GridMAT-MD. *shape* and
        *delta* are optional. The *shape* of the array is parsed from
        the filename if not provided. The spacing is set to (1,1) if
        not provided.
        """
        self.filename = filename
        self.shape = shape
        if shape is None:
            self.shape = self.parse_filename(filename)
        if delta is None:
            self.delta = tuple([1.0] * len(self.shape))   # arbitrary spacing
        else:
            self.delta = tuple(delta)

        #: 2D data from file is made available as a numpy array.
        self.array = numpy.loadtxt(self.filename).reshape(self.shape)

        #: reconstructed bins  (always start at 0,0 as the offset is lost)
        self.bins = [numpy.linspace(0,n*dx,n+1) for n,dx in zip(self.shape, self.delta)]
        #: reconstructed midpoints  (always start at 0,0 as the offset is lost)
        self.midpoints = [self._midpoints(x) for x in self.bins]

    def _midpoints(self, x):
        _x = numpy.asarray(x)
        return 0.5*(x[1:] + x[:-1])

    def parse_filename(self, filename):
        """Get dimensions from filename"""
        m = self.DATANAME.match(filename)
        if m is None:
            raise ValueError("filename %s does not appear to be a standard GridMat-MD output filename")
        return map(int, m.group('NX', 'NY'))
            
    def imshow(self, **kwargs):
        """Display data as a 2D image using :func:`pylab.imshow`."""
        import pylab
        # extent: [ None | scalars (left, right, bottom, top) ]
        extent = numpy.concatenate([self.bins[0][[0,-1]], self.bins[1][[0,-1]]])
        kwargs.setdefault('extent', extent)
        pylab.imshow(self.array, **kwargs)


    def __add__(self, *other):
        # XXX: averages would be muc nicer if we could just sum
        # XXX: these objects
        raise NotImplementedError

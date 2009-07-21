# $Id$
"""
:mod:`gromacs.analysis.plugins.gridmatmd` --- Lipid bilayer analysis helper
===========================================================================

This helper module contains code to drive ``GridMAT-MD.pl``, available
from the `GridMAT-MD`_ home page and written by WJ Allen et al
[Allen2009]_


You will need a little patch to the original version of the perl code
that is available from the GromacsWrapper distribution. Patch, rename
the script to ``GridMAT-MDx.pl``, and make sure that it can be found
on your :envvar:`PATH`.



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
    """Analysis of lipid bilayers with GridMAT-MDx."""
    
    # regular expression to pull useful information from output
    pNX = re.compile("There are (?P<NX>\d+) grid points in the X direction, "
                     "spaced every (?P<DX>[0-9.]+) nanometers")
    pNY = re.compile("There are (?P<NY>\d+) grid points in the Y direction, "
                     "spaced every (?P<DY>[0-9.]+) nanometers")
    pTOP = re.compile('The top leaflet "thickness" will be printed to (?P<FILENAME>.*_top_.*\.dat)')
    pAVG = re.compile('The average bilayer "thickness" will be printed to (?P<FILENAME>.*_average_.*\.dat)')
    pBOT = re.compile('The bottom leaflet "thickness" will be printed to (?P<FILENAME>.*_bottom_.*\.dat)')


    def __init__(self, config, gro_files):
        self.config = config
        self.filenames = gro_files   # list of filenames
        self.thickness = {}
        self.nx = self.ny = None
        self.dx = self.dy = None

        self.GridMatMD = gromacs.tools.GridMAT_MDx(self.config)

    def run_frame(self, frame):
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

        d = {}
        for name, filename in thickness.items():
            d[name] = GridMatData(filename, shape=(self.nx, self.ny),
                                  delta=(self.dx, self.dy))
        self.d = d



DATANAME = re.compile("""(?P<NX>\d+)     # number of grid points in X
                         x
                         (?P<NY>\d+)     # number of grid points in Y
                         _(top|bottom|average)_.*\.dat  # all the rest
                      """, re.VERBOSE)

class GridMatData(object):
    """Represent GridMatMD data file.

    The loaded array data is accessible as a numpy array in
    :attr:`~GridMatData.array` and bins and midpoints as
    :attr:`~GridMatData.bins` and :attr:`~GridMatData.midpoints`
    respectively.
    """
    
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
        """Get dimensions from filename, 20x20_.*.dat"""
        m = DATANAME.match(filename)
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

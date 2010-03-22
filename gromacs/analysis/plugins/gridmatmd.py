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
   :inherited-members:

.. autoclass:: Grid2D
   :members: imshow, _midpoints, __add__, __sub__, __mul__, __div__

"""
from __future__ import with_statement

import sys
import subprocess
import re
import glob
import cPickle
import numpy

import gromacs.tools

class GridMatMD(object):    
    """Analysis of lipid bilayers with GridMAT-MD.

    It requires a configuration file and a list of structure files
    (gro or pdb) as input. See the documentation (pdf_) for the format of the
    config file. Note that the *bilayer* keyword will be ignored in
    the config file.

    .. _pdf: http://www.bevanlab.biochem.vt.edu/GridMAT-MD/doc/GridMAT-MD_ug_v1.0.2.pdf
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
        """Set up GridMAT-MD analysis.

        :Arguments:
          *config* : filename
             input file for GridMAT-MD (see docs)
          *filenames* : list or glob-pattern
             list of gro or pdb files, or a glob pattern that creates
             such a list
        """
        self.config = config
        if type(filenames) is str:
            # use it as a glob pattern
            _filenames = glob.glob(filenames)
        else:
            _filenames = filenames
        #: list of filenames (gro or pdb files)        
        self.filenames = _filenames
        #: Number of bins in the *x* direction. 
        self.nx = None
        #: Numberof bins in the *y* direction.         
        self.ny = None
        #: Bin width in *x*.
        self.dx = None
        #: Bin width in *y*.
        self.dy = None
        #: Thickness results (dict of :class:`GridMatData` instances);
        #: used to store the results from the last processed frame.
        self.thickness = {}
        #: Final result: averaged thickness profiles, indexed by "top", "bottom", "average"
        self.averages = {}

        #: Instance of :class:`gromacs.tools.GridMAT_MD` which always
        #: uses *config* for input.
        self.GridMatMD = gromacs.tools.GridMAT_MD(self.config)

    def run_frame(self, frame):
        """Run GridMAT-MD on a single *frame* and return results.

        :Arguments: *frame* is a filename (gro or pdb)
        :Returns:   a dict of :class:`GridMapData` objects; the keys
                    are "top", "bottom", "average" 
        """
        
        thickness = {}
        print >> sys.stdout, "Analyzing frame %r (takes a few seconds)..." % frame
        sys.stdout.flush()
        rc,output,stderr = self.GridMatMD(frame, stdout=False)
        for line in output.split('\n'):
            print >> sys.stdout, '>>> '+line
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
        self.thickness = d  # use thickness as common block... useful for interactive work
        return d

    def run(self):
        """Run analysis on all files and average results."""

        results = {}
        for f in self.filenames:
            d = self.run_frame(f)        # <-- run GridMAT-MD
            for name,grid2d in d.items():
                try:
                    results[name].append(grid2d)
                except KeyError:
                    results[name] = [grid2d]
        
        # averages arrays and bins
        self.averages = dict([(name,numpy.mean(v)) for name,v in results.items()])

    def imshow(self, name, **kwargs):
        """Display array *name* with ``pylab.imshow``."""
        from pylab import xlabel, ylabel
        try:
            self.averages[name].imshow(**kwargs)
        except KeyError:
            raise KeyError("name must be one of %r" % self.arrays.keys())
        xlabel(r'$x/$nm')
        ylabel(r'$y/$nm')

    def save(self, name):
        """Save object as pickle."""
        with open(name, 'wb') as f:
            cPickle.dump(self, f, protocol=cPickle.HIGHEST_PROTOCOL)
        
class Grid2D(object):
    """Represents a 2D array with bin sizes.

    Addition and subtraction of grids is defined for the arrays and
    the bins.  Multiplication and division with scalars is also
    defined. Each operation returns a new :class:`Grid2D` object.

    (Actually, it should work for arrays of any dimension, not just 2D.)
    """

    def __init__(self, data, bins):
        """Initialize the Grid2D instance.

        :Arguments:
          *data*
             array data, e.g. a list of array
          *bins*
             tuple of lists of bin **edges**, one for each dimension
        """
        self.array = numpy.asarray(data)  # numpy array
        self.shape = self.array.shape
        self.ndim = len(self.shape)
        self.bins = bins
        if len(self.bins) != self.ndim:
            raise ValueError("bins must be a list of bin edges, one for each dimension, D=%(ndim)d" %
                             vars(self))
        self.midpoints = [self._midpoints(x) for x in self.bins]

    def _midpoints(self, x):
        _x = numpy.asarray(x)
        return 0.5*(_x[1:] + _x[:-1])

    def imshow(self, **kwargs):
        """Display data as a 2D image using :func:`pylab.imshow`."""
        import pylab
        # extent: [ None | scalars (left, right, bottom, top) ]
        extent = numpy.concatenate([self.bins[0][[0,-1]], self.bins[1][[0,-1]]])
        kwargs.setdefault('extent', extent)
        kwargs['origin'] = 'lower'
        pylab.imshow(self.array, **kwargs)

    def __add__(self, other):
        """Add arrays and bins (really only makes sense when averaging)."""
        if self.array.shape != other.array.shape:
            raise TypeError("arrays are of incompatible shape")
        _bins = [self.bins[dim] + other.bins[dim] for dim in xrange(len(self.bins))]
        _array = self.array + other.array
        return Grid2D(data=_array, bins=_bins)

    def __sub__(self, other):
        """Subtract other from self (also subtracts bins... which is odd but consistent)."""
        if self.array.shape != other.array.shape:
            raise TypeError("arrays are of incompatible shape")
        _bins = [self.bins[dim] - other.bins[dim] for dim in xrange(len(self.bins))]
        _array = self.array - other.array
        return Grid2D(data=_array, bins=_bins)

    def __mul__(self, x):
        """Multiply arrays (and bins) by a scalar *x*."""
        _bins = [self.bins[dim] * x for dim in xrange(len(self.bins))]
        _array = self.array * x
        return Grid2D(data=_array, bins=_bins)

    def __div__(self, x):
        """Divide arrays (and bins) by a scalar *x*."""
        _bins = [self.bins[dim] / x for dim in xrange(len(self.bins))]
        _array = self.array / x
        return Grid2D(data=_array, bins=_bins)
        


class GridMatData(Grid2D):
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

        :Arguments:
          *filename*
             2D grid as written by GridMAT-MD
          *shape*
             Shape tuple (NX, NY) of the array in filename.
          *delta*
             Tuple of bin sizes of grid (DX, DY).
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
        _array = numpy.loadtxt(self.filename).reshape(self.shape)

        #: reconstructed bins  (always start at 0,0 as the offset is lost)
        _bins = [numpy.linspace(0,n*dx,n+1) for n,dx in zip(self.shape, self.delta)]

        Grid2D.__init__(self, data=_array, bins=_bins)

        # reconstructed midpoints  (always start at 0,0 as the offset
        # is lost)        

    def parse_filename(self, filename):
        """Get dimensions from filename"""
        m = self.DATANAME.match(filename)
        if m is None:
            raise ValueError("filename %s does not appear to be a standard GridMat-MD output filename")
        return map(int, m.group('NX', 'NY'))
            

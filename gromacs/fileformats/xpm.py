# GromacsWrapper: xpm.py
# Copyright (c) 2012 Oliver Beckstein <orbeckst@gmail.com>
# Copyright (c) 2010 Tsjerk Wassenaar <tsjerkw@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.
"""
Gromacs XPM file format
=======================

Gromacs stores matrix data in the xpm file format. This implementation
of a Python reader is based on Tsjerk Wassenaar's post to gmx-users
`numerical matrix from xpm file`_ (Mon Oct 4 13:05:26 CEST 2010). This
version returns a NumPy array and can guess an appropriate dtype for
the array.

.. _numerical matrix from xpm file:
   http://lists.gromacs.org/pipermail/gmx-users/2010-October/054557.html

Classes
-------

.. autoclass:: XPM
   :members:

   .. attribute:: xvalues

      Values of on the x-axis, extracted from the xpm file.

   .. attribute:: yvalues

      Values of on the y-axis, extracted from the xpm file. These are
      in the same order as the rows in the xpm matrix. If *reverse* =
      ``False`` then this is typically a *descending* list of numbers
      (highest to lowest residue number, index number, etc). For
      *reverse* = ``True`` it is resorted accordingly.



Example: Analysing H-bonds
--------------------------

Run :func:`gromacs.g_hbond` to produce the existence map (and the log
file for the atoms involved in the bonds; the ndx file is also
useful)::

  gromacs.g_hbond(s=TPR, f=XTC, g="hbond.log", hbm="hb.xpm", hbn="hb.ndx")

Load the XPM::

  hb = XPM("hb.xpm", reverse=True)

Calculate the fraction of time that each H-bond existed::

  hb_fraction = hb.array.mean(axis=0)

Get the descriptions of the bonds::

  desc = [line.strip() for line in open("hbond.log") if not line.startswith('#')]

.. Note::

   It is important that ``reverse=True`` is set so that the rows in
   the xpm matrix are brought in the same order as the H-bond
   labels.

Show the results::

  print "\\n".join(["%-40s %4.1f%%" % p for p in zip(desc, 100*hb_fraction)])

"""

from __future__ import absolute_import, with_statement

from six.moves import range

import os, errno
import re
import warnings

import numpy

from ..exceptions import ParseError, AutoCorrectionWarning
from .. import utilities
from .convert import Autoconverter

import logging

class XPM(utilities.FileUtils):
    """Class to make a Gromacs XPM matrix available as a NumPy :class:`numpy.ndarray`.

    The data is available in the attribute :attr:`XPM.array`.

    .. Note::

       By default, the rows (2nd dimension) in the :attr:`XPM.array`
       are re-ordered so that row 0 (i.e. ``array[:,0]`` corresponds
       to the first residue/hydrogen bond/etc. The original xpm matrix
       is obtained for *reverse* = ``False``. The :class:`XPM` reader
       always reorders the :attr:`XPM.yvalues` (obtained from the xpm
       file) to match the order of the rows.

    """
    default_extension = "xpm"
    logger = logging.getLogger('gromacs.formats.XPM')
    #: compiled regular expression to parse the colors in the xpm file::
    #:
    #:   static char *gromacs_xpm[] = {
    #:   "14327 9   2 1",
    #:   "   c #FFFFFF " /* "None" */,
    #:   "o  c #FF0000 " /* "Present" */,
    #:
    #: Matches are named "symbol", "color" (hex string), and "value". "value"
    #: is typically autoconverted to appropriate values with
    #: :class:`gromacs.fileformats.convert.Autoconverter`.
    #: The symbol is matched as a `printable ASCII character`_ in the range
    #: 0x20 (space) to 0x7E (~).
    #:
    #: .. _`printable ASCII character`: http://www.danshort.com/ASCIImap/indexhex.htm
    COLOUR = re.compile("""\
            ^.*"                   # start with quotation mark
            (?P<symbol>[\x20-\x7E])# printable ASCII symbol used in the actual pixmap: 'space' to '~'
            \s+                    # white-space separated
            c\s+                   # 'c' to prefix colour??
            (?P<color>\#[0-9A-F]+) # colour as hex string (always??)
            \s*"                   # close with quotes
            \s*/\*\s*              # white space then opening C-comment /*
            "                      # start new string
            (?P<value>.*)          # description/value as free form string
            "                      # ... terminated by quotes
            """, re.VERBOSE)

    def __init__(self, filename=None, **kwargs):
        """Initialize xpm structure.

        :Arguments:
          *filename*
              read from xpm file directly
          *autoconvert*
              try to guess the type of the output array from the
              colour legend [``True``]
          *reverse*
              reverse rows (2nd dimension): re-orders the rows so that
              the first row corresponds e.g. to the first residue or
              first H-bonds and not the last) [``True``]
        """
        self.autoconvert = kwargs.pop("autoconvert", True)
        self.reverse = kwargs.pop("reverse", True)
        self.__array = None
        super(XPM, self).__init__(**kwargs)  # can use kwargs to set dict! (but no sanity checks!)

        if filename is not None:
            self._init_filename(filename)
            self.read(filename)


    def to_df(self):
        import pandas as pd

        # Add Time to the data as column
        data = numpy.vstack((self.xvalues, self.array.T)).T

        # Column names are resids
        df = pd.DataFrame(data, columns=["Time"]+ list(self.yvalues))

        # Converts Time to a numeric type
        df['Time'] = pd.to_numeric(df['Time'])
        return df

    @property
    def array(self):
        """XPM matrix as a :class:`numpy.ndarray`.

        The attribute itself cannot be assigned a different array but
        the contents of the array can be modified.
        """
        return self.__array

    def read(self, filename=None):
        """Read and parse mdp file *filename*."""
        self._init_filename(filename)
        self.parse()

    def parse(self):
        """Parse the xpm file and populate :attr:`XPM.array`."""
        with utilities.openany(self.real_filename) as xpm:
            # Read in lines until we find the start of the array
            meta = [xpm.readline()]
            while not meta[-1].startswith("static char *gromacs_xpm[]"):
                meta.append(xpm.readline())

            # The next line will contain the dimensions of the array
            dim = xpm.readline()
            # There are four integers surrounded by quotes
            # nx: points along x, ny: points along y, nc: ?, nb: stride x
            nx, ny, nc, nb = [int(i) for i in self.unquote(dim).split()]

            # The next dim[2] lines contain the color definitions
            # Each pixel is encoded by dim[3] bytes, and a comment
            # at the end of the line contains the corresponding value
            colors = dict([self.col(xpm.readline()) for i in range(nc)])


            if self.autoconvert:
                autoconverter = Autoconverter(mode="singlet")
                for symbol, value in colors.items():
                    colors[symbol] = autoconverter.convert(value)
                self.logger.debug("Autoconverted colours: %r", colors)

            # make an array containing all possible values and let numpy figure out the dtype
            dtype = numpy.array(colors.values()).dtype
            self.logger.debug("Guessed array type: %s", dtype.name)

            # pre-allocate array
            data = numpy.zeros((int(nx/nb), ny), dtype=dtype)

            self.logger.debug("dimensions: NX=%d NY=%d strideX=%d (NC=%d) --> (%d, %d)",
                              nx, ny, nb, nc, nx/nb, ny)

            iy = 0
            xval = []
            yval = []
            autoconverter = Autoconverter(mode="singlet")
            for line in xpm:
                if line.startswith("/*"):
                    # lines '/* x-axis:' ... and '/* y-axis:' contain the
                    # values of x and y coordinates
                    s = self.uncomment(line).strip()
                    if s.startswith('x-axis:'):
                        xval.extend([autoconverter.convert(x) for x in s[7:].split()])
                    elif s.startswith('y-axis:'):
                        yval.extend([autoconverter.convert(y) for y in s[7:].split()])
                    continue
                s = self.unquote(line)
                # Joao M. Damas <jmdamas@itqb.unl.pt> suggests on gmx-users (24 Oct 2014)
                # that the next line should read:
                #
                #  data[:, iy]  =  [colors[j[k:k+nb]] for k in range(0,nx*nb,nb)]
                #
                # "if one is using higher -nlevels for the .xpm construction (in g_rms, for example)"
                # However, without a test case I am not eager to change it right away so in
                # case some faulty behavior is discovered with the XPM reader then this comment
                # might be helpful. --- Oliver 2014-10-25
                data[:, iy] = [colors[s[k:k+nb]] for k in range(0,nx,nb)]
                self.logger.debug("read row %d with %d columns: '%s....%s'",
                                  iy, data.shape[0], s[:4], s[-4:])
                iy += 1  # for next row

        self.xvalues = numpy.array(xval)
        if self.reverse:
            self.logger.debug("reversed row order, reverse=%r", self.reverse)
            self.__array = data[:, ::-1]
            self.yvalues = numpy.array(yval)
        else:
            self.__array = data
            self.yvalues = numpy.array(yval)[::-1]  # must reverse y-values to match!

    @staticmethod
    def unquote(s):
        """Return string *s* with quotes ``"`` removed."""
        return s[1+s.find('"'):s.rfind('"')]

    @staticmethod
    def uncomment(s):
        """Return string *s* with C-style comments ``/*`` ... ``*/`` removed."""
        return s[2+s.find('/*'):s.rfind('*/')]


    def col(self, c):
        """Parse colour specification"""
        m = self.COLOUR.search(c)
        if not m:
            self.logger.fatal("Cannot parse colour specification %r.", c)
            raise ParseError("XPM reader: Cannot parse colour specification {0!r}.".format(c))
        value = m.group('value')
        color = m.group('symbol')
        self.logger.debug("%s: %s %s\n", c.strip(), color, value)
        return color, value


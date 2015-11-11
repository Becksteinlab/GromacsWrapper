# GromacsWrapper: formats.py
# Copyright (c) 2009-2011 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.
"""
GRO structure file format
=========================

**Not implemented.**

..SeeAlso:: MDAnalysis_ has a working Python `GRO reader`_.

.. _MDAnalysis: http://mdanalysis.googlecode.com
.. _GRO reader: http://mdanalysis.googlecode.com/svn/trunk/doc/html/documentation_pages/coordinates/GRO.html

.. autoclass:: GRO
   :members:

"""


from __future__ import absolute_import, with_statement

import os, errno
import re
import warnings

import numpy

from ..exceptions import ParseError, AutoCorrectionWarning
from .. import utilities
try:
    from collections import OrderedDict as odict
except ImportError:
    from .odict import odict

import logging

class GRO(utilities.FileUtils):
    """Class that represents a GROMOS (gro) structure file.


    File format:
    """
    default_extension = "gro"
    logger = logging.getLogger('gromacs.formats.GRO')

    def __init__(self, **kwargs):

        raise NotImplementedError

        filename = kwargs.pop('filename',None)
        super(GRO, self).__init__(**kwargs)

        if not filename is None:
            self._init_filename(filename)
            self.read(filename)

    def read(self, filename=None):
        """Read and parse index file *filename*."""
        self._init_filename(filename)

        with open(self.real_filename) as gro:
            pass



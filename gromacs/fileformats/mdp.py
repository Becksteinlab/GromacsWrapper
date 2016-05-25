# GromacsWrapper: formats.py
# Copyright (c) 2009-2011 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.
"""
Gromacs parameter MDP file format
=================================

The `.mdp file`_ contains a list of keywords that are used to set up a
simulation with :class:`~gromacs.tools.Grompp`. The class :class:`MDP`
parses this file and provides access to the keys and values as ordered
dictionary.

.. _`.mdp file`: http://www.gromacs.org/Documentation/File_Formats/.mdp_File

.. autoclass:: MDP
   :members:
"""


from __future__ import absolute_import, with_statement

import os, errno
import re
import warnings

import numpy


from ..exceptions import ParseError, AutoCorrectionWarning
from .. import utilities
from collections import OrderedDict as odict

import logging

class MDP(odict, utilities.FileUtils):
    """Class that represents a Gromacs mdp run input file.

    The MDP instance is an ordered dictionary.

      - *Parameter names* are keys in the dictionary.
      - *Comments* are sequentially numbered with keys Comment0001,
        Comment0002, ...
      - *Empty lines* are similarly preserved as Blank0001, ....

    When writing, the dictionary is dumped in the recorded order to a
    file. Inserting keys at a specific position is not possible.

    Currently, comments after a parameter on the same line are
    discarded. Leading and trailing spaces are always stripped.

    .. SeeAlso:: For editing a mdp file one can also use
                :func:`gromacs.cbook.edit_mdp` (which works like a
                poor replacement for sed).
    """
    default_extension = "mdp"
    logger = logging.getLogger('gromacs.formats.MDP')

    COMMENT = re.compile("""\s*;\s*(?P<value>.*)""")   # eat initial ws
    # see regex in cbook.edit_mdp()
    PARAMETER = re.compile("""
                            \s*(?P<parameter>[^=]+?)\s*=\s*  # parameter (ws-stripped), before '='
                            (?P<value>[^;]*)                # value (stop before comment=;)
                            (?P<comment>\s*;.*)?            # optional comment
                            """, re.VERBOSE)

    def __init__(self, filename=None, autoconvert=True, **kwargs):
        """Initialize mdp structure.

        :Arguments:
          *filename*
              read from mdp file
          *autoconvert* : boolean
              ``True`` converts numerical values to python numerical types;
              ``False`` keeps everything as strings [``True``]
          *kwargs*
              Populate the MDP with key=value pairs. (NO SANITY CHECKS; and also
              does not work for keys that are not legal python variable names such
              as anything that includes a minus '-' sign or starts with a number).
        """
        super(MDP, self).__init__(**kwargs)  # can use kwargs to set dict! (but no sanity checks!)

        self.autoconvert = autoconvert

        if not filename is None:
            self._init_filename(filename)
            self.read(filename)

    def _transform(self, value):
        if self.autoconvert:
            return utilities.autoconvert(value)
        else:
            return value

    def read(self, filename=None):
        """Read and parse mdp file *filename*."""
        self._init_filename(filename)

        def BLANK(i):
            return "B%04d" % i
        def COMMENT(i):
            return "C%04d" % i

        data = odict()
        iblank = icomment = 0
        with open(self.real_filename) as mdp:
            for line in mdp:
                line = line.strip()
                if len(line) == 0:
                    iblank += 1
                    data[BLANK(iblank)] = ''
                    continue
                m = self.COMMENT.match(line)
                if m:
                    icomment += 1
                    data[COMMENT(icomment)] = m.group('value')
                    continue
                # parameter
                m = self.PARAMETER.match(line)
                if m:
                    # check for comments after parameter?? -- currently discarded
                    parameter = m.group('parameter')
                    value =  self._transform(m.group('value'))
                    data[parameter] = value
                else:
                    errmsg = '%(filename)r: unknown line in mdp file, %(line)r' % vars()
                    self.logger.error(errmsg)
                    raise ParseError(errmsg)

        super(MDP,self).update(data)


    def write(self, filename=None, skipempty=False):
        """Write mdp file to *filename*.

        :Keywords:
           *filename*
               output mdp file; default is the filename the mdp
               was read from
           *skipempty* : boolean
               ``True`` removes any parameter lines from output that
               contain empty values [``False``]

        .. Note:: Overwrites the file that the mdp was read from if no
                  *filename* supplied.
        """

        with open(self.filename(filename, ext='mdp'), 'w') as mdp:
            for k,v in self.items():
                if k[0] == 'B':        # blank line
                    mdp.write("\n")
                elif k[0] == 'C':      # comment
                    mdp.write("; %(v)s\n" % vars())
                else:                  # parameter = value
                    if skipempty and (v == '' or v is None):
                        continue
                    mdp.write("%(k)s = %(v)s\n" % vars())


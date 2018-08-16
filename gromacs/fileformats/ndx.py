# GromacsWrapper: formats.py
# Copyright (c) 2009-2011 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
Gromacs NDX index file format
=============================

The `.ndx file`_ contains lists of atom indices that are grouped in
sections by *group names*. The classes :class:`NDX` and
:class:`uniqueNDX` can parse such ndx files and provide convenient
access to the individual groups.

.. _`.ndx file`: http://www.gromacs.org/Documentation/File_Formats/.ndx_File

.. autoclass:: NDX
   :members:

.. autoclass:: uniqueNDX
   :members:

.. autoclass:: IndexSet
"""

from __future__ import absolute_import, with_statement

from six.moves import range

import os, errno
import re
import warnings
import operator

import numpy

from ..exceptions import ParseError, AutoCorrectionWarning
from .. import utilities
from collections import OrderedDict as odict

import logging

class NDX(odict, utilities.FileUtils):
    """Gromacs index file.

    Represented as a ordered dict where the keys are index group names and
    values are numpy arrays of atom numbers.

    Use the :meth:`NDX.read` and :meth:`NDX.write` methods for
    I/O. Access groups by name via the :meth:`NDX.get` and
    :meth:`NDX.set` methods.

    Alternatively, simply treat the :class:`NDX` instance as a
    dictionary. Setting a key automatically transforms the new value
    into a integer 1D numpy array (*not* a set, as would be the
    :program:`make_ndx` behaviour).

    .. Note::

       The index entries themselves are ordered and can contain
       duplicates so that output from NDX can be easily used for
       :program:`g_dih` and friends. If you need set-like behaviour
       you will have do use :class:`gromacs.formats.uniqueNDX` or
       :class:`gromacs.cbook.IndexBuilder` (which uses
       :program:`make_ndx` throughout).

    **Example**

      Read index file, make new group and write to disk::

        ndx = NDX()
        ndx.read('system.ndx')
        print ndx['Protein']
        ndx['my_group'] = [2, 4, 1, 5]   # add new group
        ndx.write('new.ndx')

      Or quicker (replacing the input file ``system.ndx``)::

        ndx = NDX('system')          # suffix .ndx is automatically added
        ndx['chi1'] = [2, 7, 8, 10]
        ndx.write()

    """
    default_extension = "ndx"

    # match:  [ index_groupname ]
    SECTION = re.compile("""\s*\[\s*(?P<name>\S.*\S)\s*\]\s*""")

    #: standard ndx file format: 15 columns
    ncol = 15
    #: standard ndx file format: '%6d'
    format = '%6d'

    def __init__(self, filename=None, **kwargs):
        super(NDX, self).__init__(**kwargs)  # can use kwargs to set dict! (but no sanity checks!)

        if filename is not None:
            self._init_filename(filename)
            self.read(filename)

    def read(self, filename=None):
        """Read and parse index file *filename*."""
        self._init_filename(filename)

        data = odict()
        with open(self.real_filename) as ndx:
            current_section = None
            for line in ndx:
                line = line.strip()
                if len(line) == 0:
                    continue
                m = self.SECTION.match(line)
                if m:
                    current_section = m.group('name')
                    data[current_section] = []  # can fail if name not legal python key
                    continue
                if current_section is not None:
                    data[current_section].extend(map(int, line.split()))

        super(NDX,self).update(odict([(name, self._transform(atomnumbers))
                                     for name, atomnumbers in data.items()]))

    def write(self, filename=None, ncol=ncol, format=format):
        """Write index file to *filename* (or overwrite the file that the index was read from)"""
        with open(self.filename(filename, ext='ndx'), 'w') as ndx:
            for name in self:
                atomnumbers = self._getarray(name)  # allows overriding
                ndx.write('[ {0!s} ]\n'.format(name))
                for k in range(0, len(atomnumbers), ncol):
                    line = atomnumbers[k:k+ncol].astype(int)   # nice formatting in ncol-blocks
                    n = len(line)
                    ndx.write((" ".join(n*[format])+'\n') % tuple(line))
                ndx.write('\n')

    def get(self, name):
        """Return index array for index group *name*."""
        return self[name]

    def set(self, name, value):
        """Set or add group *name* as a 1D numpy array."""
        self[name] = value

    def size(self, name):
        """Return number of entries for group *name*."""
        return len(self[name])

    @property
    def groups(self):
        """Return a list of all groups."""
        return self.keys()

    @property
    def sizes(self):
        """Return a dict with group names and number of entries,"""
        return {name: len(atomnumbers) for name, atomnumbers in self.items()}

    @property
    def ndxlist(self):
        """Return a list of groups in the same format as  :func:`gromacs.cbook.get_ndx_groups`.

        Format:
           [ {'name': group_name, 'natoms': number_atoms, 'nr':  # group_number}, ....]
        """
        return [{'name': name, 'natoms': len(atomnumbers), 'nr': nr+1} for
                nr,(name,atomnumbers) in enumerate(self.items())]

    def _getarray(self, name):
        """Helper getter that is used in write().
        Override when using a _transform that stores something that
        cannot be indexed, e.g. when using set()s.
        """
        return self[name]

    def _transform(self, v):
        """Transform input to the stored representation.

        Override eg with ``return set(v)`` for index lists as sets.
        """
        return numpy.ravel(v).astype(int)

    def __setitem__(self, k, v):
        super(NDX, self).__setitem__(k, self._transform(v))

    def setdefault(*args,**kwargs):
        raise NotImplementedError


class IndexSet(set):
    """set which defines '+' as union (OR) and '-' as intersection  (AND)."""
    def __add__(self, x):
        return self.union(x)
    def __sub__(self, x):
        return self.intersection(x)


class uniqueNDX(NDX):
    """Index that behaves like make_ndx, i.e. entries behaves as sets,
    not lists.

    The index lists behave like sets:
    - adding sets with '+' is equivalent to a logical OR: x + y == "x | y"
    - subtraction '-' is AND: x - y == "x & y"
    - see :meth:`~gromacs.formats.join` for ORing multiple groups (x+y+z+...)

    **Example** ::

       I = uniqueNDX('system.ndx')
       I['SOLVENT'] = I['SOL'] + I['NA+'] + I['CL-']

    """

    def join(self, *groupnames):
        """Return an index group that contains atoms from all  *groupnames*.

        The method will silently ignore any groups that are not in the
        index.

        **Example**

        Always make a solvent group from water and ions, even if not
        all ions are present in all simulations::

           I['SOLVENT'] = I.join('SOL', 'NA+', 'K+', 'CL-')
        """
        return self._sum([self[k] for k in groupnames if k in self])

    def _sum(self, sequence):
        return reduce(operator.add, sequence)

    def _transform(self, v):
        return IndexSet(v)

    def _getarray(self, k):
        return numpy.sort(numpy.fromiter(self[k],dtype=int,count=len(self[k])))



# or use list of these?
# class IndexGroup(dict):
#     def __init__(self, groupnumber=None, name="empty", atomnumbers=None, **kwargs):
#         atomnumbers = atomnumbers or []
#         _atomnumbers = numpy.asarray(atomnumbers).astype(int)
#         super(IndexGroup, self).__init__(name=str(name),
#                                          atomnumbers=_atomnumbers,
#                                          nr=groupnumber)

# GromacsWrapper: formats.py
# Copyright (c) 2009-2011 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.
"""
Gromacs topology file (ITP) parser
==================================

.. versionadded:: 0.2.5

Basic reading, manipulating and writing works but can still be a bit
awkward. See documentation for :class:`ITP`.

.. warning:: The code is still *experimental* and the interface
             will probably change somewhat until release 0.3.

.. rubric:: Limitations

- merges multiple dihedral sections into one
- probably fails for FEP itps (TODO)
- does not reproduce positions of comment-only lines in sections
  (in fact, currently we do not even write them out again, only trailing
  line comments are kept and the header before the first section)
- fails with preprocessor directives such as ``#ifdef POSRES``

User classes
------------

.. autoclass:: ITP
   :members:

Developer reference
-------------------

Class :class:`ITP` makes use of the following classes to parse and
represent a topology as stored in an itp file.

Base classes and utility functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: ITPsection
   :members:

.. autoclass:: ITPdata
   :members:

.. autoclass:: OneLineBuffer
   :members:

.. autofunction:: flatten

Classes corresponding to ITP sections
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ITP file is parsed with a simple recursive descending depth-first
algorithm. Classes contain a dictionary of sub-section parsers.

.. autoclass:: Header
   :members:

.. autoclass:: Moleculetype
   :members:

.. autoclass:: Atoms
   :members:

.. autoclass:: Bonds
   :members:

.. autoclass:: Angles
   :members:

.. autoclass:: Dihedrals
   :members:

.. autoclass:: Pairs
   :members:

.. autoclass:: Dummy
   :members:

"""

from __future__ import with_statement
import os, errno
import re
import warnings

import numpy

from gromacs import ParseError, AutoCorrectionWarning
import gromacs.utilities as utilities
from gromacs.odict import odict

import logging



class ITPsection(object):
    """Class to parse and store ITP data.

    - :attr:`data` contains the data of the section
    - :attr:`sections` holds the sub-parsers (including their own :attr:`data`;
      this :class:`odict` mirrors the parse tree.

      Note that sections at the same level are merged.

    - :attr:`parsers` contains a dict of the all the parser classes
      (keyed by section name) that are included in *this* section;
      typically overriden.

     .. versionadded:: 0.2.5
    """
    COMMENT = re.compile("""\s*;\s*(?P<comment>.*)""")       # eat initial ws
    SECTION = re.compile("""^\s*\[\s*(?P<name>\S+)\s*\]""")  # [ name ]
    DATA = re.compile("""\
                      ^\s*                    # ignore initial ws
                      (?P<data>[^;]+)         # whole line except trailing comment
                      (;\s*                   # eat leading white space
                        (?P<comment>.*)       # possibly a comment
                      )?
                      """, re.VERBOSE)
    name = None

    def __init__(self, parent, **kwargs):
        self.parent = parent        # parent section
        self.logger = parent.logger
        self.sections = odict()
        self.__data = kwargs.pop("data", [])
        self.comments = kwargs.pop("comments", [])
        self.parsers = {}

    @property
    def data(self):
        return self.__data
    def section(self):
        # incomplete, override!
        return "[ %s ]" % self.name
    def parse(self, line):
        pass   # implement in subclasses
    def __len__(self):
        try:
            return len(self.data)
        except (TypeError, ValueError):
            return 0
    def __str__(self):
        """Print full section as seen in an ITP file"""
        return self.section()
    def __repr__(self):
        return "<ITP::%s (%d entries)>" % (self.__class__.__name__, len(self))

    def parse(self, stream):
        current_section = self.name
        while True:
            try:
                line = stream.next().strip()
            except (StopIteration, EOFError):
                break                  # done reading the input

            if len(line) == 0:
                continue               # skip empty lines for now

            m = self.SECTION.match(line)
            if m:
                if current_section == m.group('name'):
                    self.logger.info("Merging [ %(current_section)s ] sections!", vars())
                    # Note: If we store sections in a list instead of an odict then
                    #       we can keep them separate and don't have to merge here.
                    continue
                # switch to new section and parser
                current_section = m.group('name')
                try:
                    parser = self.parsers[current_section](self)
                except KeyError:
                    # no parsers in this section...  but need to
                    # unread the line so that the next parser can get
                    # the section...
                    stream.unread()
                    break
                self.sections[current_section] = parser
                self.logger.debug("Parsing section %(current_section)r", vars())

                parser.parse(stream)
            else:
                self.process(line)

    def walk_sections(self, func=None, flat=True):
        """Generator that applies *func* to each section and subsections.

        Traverses all sections depth first. Applies

        .. function:: func(name, section)

           User supplied function that uses both the string *name* and
           the :class:`ITPsection` instance and returns a result. Note
           that any lists and lists of lists will be flattened so it
           mostly makes sense to return a string or similar.

        at each level and iterates over the individual results. With
        *flat* = ``True`` the resulting list is returned flat
        (default), otherwise the list represents the tree structure of
        the itp file.

        By default, *func* is ``str(section)``.
        """
        if func is None:
            def func(name, section):
                return str(section)
        values = [func(self.name, self)]
        for name, section in self.sections.items():
            values.append(section.walk_sections(func, flat=flat))
        if flat:
            return flatten(values)
        else:
            return values

class ITPdata(ITPsection):
    """Class to represent most ITP sections.

    1) parse line by line and store data and comments
    2) format data as recarray; uses up columns from left to right
    3) manipulate data in the recarray or write section with :meth:`section`

    .. versionadded:: 0.2.5
    """
    dtypes = []  # override, define columns [("name", "type"), ...] from left to right
    fmt = []     # override, define output fmt "%6d" (used in section())
    column_comment = "" # line to be output after section header

    def __init__(self, *args, **kwargs):
        super(ITPdata, self).__init__(*args, **kwargs)
        self.__records = []   # [((atomnr, atomtype, ...), comment), ...]
        self.__data = None    # recarray

    def _create_recarray(self):
        """Build a recarray from parsed data."""
        # build empty record array and fill line by line to accomodate
        # lines with variable number of entries
        numcols = [len(record)  for record, comment in self.__records]
        nmax = numpy.max(numcols)
        dtype = self.dtypes[:nmax] + [("comment", "S128")]
        a = numpy.recarray((len(self.__records),), dtype=dtype)
        for i, (record, comment) in enumerate(self._canonical_records(nmax=nmax)):
            a[i] = record + (comment,)
        return a

    def _canonical_records(self, nmax=None):
        """Bring records to shape with nmax columns.

        Empty columns are filled with ``None``.
        """
        if nmax is None:
            nmax = numpy.max([len(record)  for record, comment in self.__records])
        records = []
        for record, comment in self.__records:
            new_col = nmax - len(record)
            records.append((tuple(record) + new_col * (None,), comment))
        return records

    @property
    def data(self):
        """atom data as a :class:`numpy.rec.array`

        The data inside the array can be changed.
        """
        if self.__data is None:
            self.__data = self._create_recarray()
        return self.__data

    def process(self, line):
        if len(line.strip()) == 0:
            return     # skip empty lines

        m = self.COMMENT.match(line)
        if m:
            self.comments.append(m.group('comment'))
            return

        m = self.DATA.match(line)
        if m:
            comment = m.group('comment') or ""
            record = m.group('data').split()
            self.__records.append((record, comment))
            return

        self.logger.warn("[%s] not parsing line: %r", self.name, line)

    def section(self):
        def clean_records(data):
            for rec in data:
                yield tuple([x for x in rec if x is not None])

        lines = ["[ %s ]" % self.name]
        lines.append(self.column_comment)

        for rec in clean_records(self.data):
            numcols = len(rec) - 1    # subtract 1 because comment is always last field in record
            fmt = " ".join(self.fmt[:numcols])
            if rec[-1]:
                fmt += " ; %s"
                line = fmt % tuple(rec)
            else:
                line = fmt % tuple(rec)[:-1]
            lines.append(line)

        return "\n".join(lines) + "\n"

class Header(ITPsection):
    """The customary comments section before the first section.

    .. versionadded:: 0.2.5
    """
    name = "header"

    def __init__(self, *args, **kwargs):
        super(Header, self).__init__(*args, **kwargs)
        self.parsers = {'moleculetype': Moleculetype}

    def process(self, line):
        self.data.append(line)
    def section(self):
        return "\n".join(self.data) + "\n"

class Moleculetype(ITPsection):
    """ITP ``[ moleculetype ]`` section

    ::
       [ moleculetype ]
       ; Name      nrexcl
       5FH              3

    .. versionadded:: 0.2.5
    """
    name = "moleculetype"

    def __init__(self, *args, **kwargs):
        kwargs['data'] = odict()
        super(Moleculetype, self).__init__(*args, **kwargs)
        self.parsers = {
            'atoms': Atoms,
            'bonds': Bonds,
            'angles': Angles,
            'dihedrals': Dihedrals,
            'pairs': Pairs,
            'dummy': Dummy,
            }

    def process(self, line):
        if len(line.strip()) == 0:
            return     # skip empty lines
        m = self.COMMENT.match(line)
        if m:
            self.comments.append("%(comment)s" % m.groupdict())
            return
        fields = line.split()
        try:
            self.data['name'] = fields[0]
            self.data['nrexcl'] = int(fields[1])
        except Exception, err:
            msg = "Failed to parse [moleculetype] section: line: %r\n" % line
            msg += "Exception: %r" % err
            self.logger.error(msg)
            raise

    def set_moleculename(self, molname):
        """Set the molecule name to *molname*.

        1. changes the name in ``[moleculetype]`` section
        2. changes the resname in the whole ``[atoms]`` section
        """
        self.data['name'] = molname
        self.sections['atoms'].set_resname(molname)

    def section(self):
        # currently without user comments
        return "[ %s ]\n; Name      nrexcl\n" % self.name + \
            "%(name)-10s  %(nrexcl)d" % self.data  + "\n"

    def __repr__(self):
        return "<ITP::moleculetype %(name)s nrexcl=%(nrexcl)d>" % self.data


class Atoms(ITPdata):
    """ITP ``[ atoms ]`` section

    ::
       [ atoms ]
       ; atomnr  atomtype   resnr  resname  atomname  chargegrp   charge       mass
              1  opls_145       1      5FH        C7          7   -0.115   12.01100 ; CA # Benzene C - 12 site JACS,112,4768-90. Use #145B for biphenyl
              2  opls_146       1      5FH       H71          7    0.115    1.00800 ; HA # Benzene H - 12 site.


    .. Warning:: Currently fails to correctly parse ITPs containing multiple
                 states for an atom as used in FEP calculations. The parser only
                 handles files with the same number of columns for each entry.

    The data itself are stored in :attr:`data`, a :class:`numpy.rec.array`.
    :attr:`data` is a managed attribute but the values inside can be
    changed. Fields in the array have a maximum size (in particular, comments
    are cut off after 128 characters).

    .. versionadded:: 0.2.5
    """
    name = "atoms"
    #: :class:`numpy.dtype` columns for the data
    dtypes = [("atomnr", "i4"), ("atomtype", "S10"), ("resnr", "i4"), ("resname", "S4"),
              ("atomname", "S4"), ("chargegrp", "i4"),
              ("charge", "f8"), ("mass", "f8"),
              ("atomtypeB", "S10"),
              ("chargeB", "f8"), ("massB", "f8"),
              # extra last column ("comment", "S128")
              ]
    #: output format (space separated), same ordering as :attr:`columns`
    fmt = ["%8d", "%9s", "%7d", "%8s",
           "%9s", "%10d",
           "%8.3f", "%10.5f",
           "%9s",
           "%8.3f", "%10.5f",
           ]
    column_comment = "; atomnr  atomtype   resnr  resname  atomname  chargegrp   charge       mass"
    # XXX: should get the column_comment from dtypes and fmt!

    def set_resname(self, name):
        """Changes the resname to *name* for all atom entries"""
        self.data.resname = name  # changes all resname entries at the same time :-)

class Bonds(ITPdata):
    """ITP ``[ bonds ]`` section
    ::
       [ bonds ]
       ; ai   aj  funct  r  k
         1    2      1 ; CA-HA # PHE, etc.
         1    3      1 ; CA-CA # TRP,TYR,PHE

    .. versionadded:: 0.2.5
    """
    name = "bonds"
    #: :class:`numpy.dtype` columns for the data
    dtypes = [("ai", "i4"), ("aj", "i4"), ("func", "i4"),
              ("b0", "f8"), ("kb", "f8"),
              # XXX: FEP columns?
              # extra last column ("comment", "S128")
              ]
    #: output format (space separated), same ordering as :attr:`columns`
    fmt = ["%4d", "%4d", "%6d",
           "%10.5f", "%10.1f",
           ]
    column_comment = "; ai   aj  funct  r  k"


class Angles(ITPdata):
    """ITP ``[ angles ]`` section

    ::
       [ angles ]
       ; ai   aj   ak  funct  theta   cth
          2    1    3      1 ; HA-CA-CA #
          2    1   11      1 ; HA-CA-CA #
          3    1   11      1 ; CA-CA-CA # PHE(OL)

    .. versionadded:: 0.2.5
    """
    name = "angles"
    #: :class:`numpy.dtype` columns for the data
    dtypes = [("ai", "i4"), ("aj", "i4"), ("ak", "i4"), ("func", "i4"),
              ("theta0", "f8"), ("cth", "f8"),
              # XXX: FEP columns?
              # extra last column ("comment", "S128")
              ]
    #: output format (space separated), same ordering as :attr:`columns`
    fmt = ["%4d", "%4d", "%4d", "%6d",
           "%9.3f", "%10.3f",
           ]
    column_comment = "; ai   aj   ak  funct  theta   cth"


class Dihedrals(ITPdata):
    """ITP ``[ dihedrals ]`` section

    ::
       [ dihedrals ]
       ; ai   aj   ak   al  funct   C0  ...  C5
          2    1    3    4      3 ; HA-CA-CA-HA # aromatic ring (X-CA-CA-X generic proper dihedral)
          2    1    3    5      3 ; HA-CA-CA-CA # aromatic ring (X-CA-CA-X generic proper dihedral)

    .. Note:: Multiple dihedral sections are currently merged into a single
             dihedral section. This means that if you separated proper and
             improper dihedrals into two different sections then they will now
             appear one after another under a single ``[ dihedrals ]`` section.

    .. versionadded:: 0.2.5
    """
    # NOTE: XXX this can appear multiple times (propers and impropers!)

    name = "dihedrals"
    #: :class:`numpy.dtype` columns for the data
    dtypes = [("ai", "i4"), ("aj", "i4"), ("ak", "i4"), ("al", "i4"), ("func", "i4"),
              #("c0", "f8"), ("c1", "f8"), ("c2", "f8"), ("c3", "f8"), ("c4", "f8"), ("c5", "f8"),
              ("c0", object), ("c1", object), ("c2", object), ("c3", object), ("c4", object), ("c5", object),
              # XXX: FEP columns?
              # extra last column ("comment", "S128")
              ]
    #: output format (space separated), same ordering as :attr:`columns`
    fmt = ["%4d", "%4d", "%4d", "%4d", "%5d",
           #"%f9.5", "%f9.5", "%f9.5", "%f9.5", "%f9.5", "%f9.5",
           "%9s", "%9s", "%9s", "%9s", "%9s", "%9s",
           ]
    column_comment = "; ai   aj   ak   al  funct   C0  ...  C5"


class Pairs(ITPdata):
    """ITP [ pairs ] section

    ::
       [ pairs ]
       ; ai   aj  funct
          1    6      1
          1    7      1

    .. versionadded:: 0.2.5
    """
    name = "pairs"
    dtypes = [("ai", "i4"), ("aj", "i4"), ("func", "i4"), ]
    fmt = ["%4d", "%4d", "%6d"]
    column_comment = "; ai   aj  funct"


class Dummy(ITPsection):
    """Not implemented, just ignores everything"""
    name = "NOT_IMPLEMENTED"
    def process(self, line):
        pass


class ITP(utilities.FileUtils):
    """Class that represents a Gromacs topology itp file.

    The file is represented by the sections (and subsections) of the
    attribute :attr:`ITP.sections`.

    Data of (sub)sections is stored in :attr:`data`, e.g. ::

        itp.sections['header'].sections['moleculetype'].sections['atoms'].data

    contains the list of atoms for molecule ::

        itp.sections['header'].sections['moleculetype'].data['name']

    Data attributes are typically (ordered) dictionaries or NumPy recarrays.

    The top-level section is called *"header"* and simply contains all
    comment lines before the first real parsed section.

    .. versionadded:: 0.2.5
    """
    default_extension = "itp"
    logger = logging.getLogger('gromacs.formats.ITP')

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
        super(ITP, self).__init__(**kwargs)

        self.autoconvert = autoconvert
        self.sections = odict()
        self.parsers  = {
            'header': Header,
            'dummy': Dummy,
            }

        if not filename is None:
            self._init_filename(filename)
            self.read(filename)

    def _transform(self, value):
        if self.autoconvert:
            return autoconvert(value)
        else:
            return value

    def read(self, filename=None):
        """Read and parse itp file *filename*."""
        self._init_filename(filename)

        with open(self.real_filename) as itp:
            stream = OneLineBuffer(itp.next)
            self.parse(stream)

    def parse(self, stream):
        """Simple recursive descent parsing of ITP.

        *stream* must have a ``next()`` and an ``unread()`` method. A
        call must supply the next line of input or raise
        :exc:`EOFError` or :exc:`StopIteration`.
        """
        # store by section

        current_section = "header"
        parser = self.parsers[current_section](self)  # argument is the parent!
        self.sections[current_section] = parser

        parser.parse(stream)

    def walk_sections(self, func=None, flat=True):
        """Generator that applies *func* to each section.

        Traverses all sections depth first. Applies

        .. function:: func(name, section)

           User supplied function that uses both the string *name* and
           the :class:`ITPsection` instance and returns a result. Note
           that any lists and lists of lists will be flattened so it
           mostly makes sense to return a string or similar.

        at each level and iterates over the individual results. With
        *flat* = ``True`` the resulting list is returned flat
        (default), otherwise the list represents the tree structure of
        the itp file.

        By default, *func* is ``str(section)``. Thus, the whole itp
        file can be printed out with ::

          print "\n".join(itp.walk_sections())

        (See :meth:`write` for writing out the ITP.)
        """
        # method copied *almost* verbatim from ITPsection.walk_sections()
        if func is None:
            def func(name, section):
                return str(section)
        values = []
        for name, section in self.sections.items():
            values.append(section.walk_sections(func, flat=flat))
        if flat:
            return flatten(values)
        else:
            return values

    def write(self, filename):
        """Write ITP file to *filename*"""
        with open(filename, "w") as itp:
            itp.write(str(self))

    def __str__(self):
        """Printable representation: whole ITP file."""
        return "\n".join(self.walk_sections())


class OneLineBuffer(object):
    """Wrapper around anything with a ``next()`` method, to provide an ``unread()``.

    Provide the ``next`` method (or any other function that returns
    one line of input) as the argument *getline*.

    Implements a one-line buffer, :meth:`unread` can only move the
    stream one line back (repeatedly calling :meth:`unread` does not
    do anything).

    .. versionadded:: 0.2.5
    """
    def __init__(self, getline):
        self.buffer = None
        self.getline = getline
        self._reread_last_line = False
    def next(self):
        """Return next line (or previous line if :meth:`unread` was called)"""
        if not self._reread_last_line:
            self.buffer = self.getline()
        self._reread_last_line = False
        return self.buffer
    def unread(self):
        """Reset stream to previous line, so that next call to :meth:`next` rereads"""
        self._reread_last_line = True


import collections  # python 2.6?

def flatten(l):
    """Generator for flattening a nested list *l*.

    http://stackoverflow.com/questions/2158395/flatten-an-irregular-list-of-lists-in-python

    .. versionadded:: 0.2.5
    """
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el, basestring):
            for sub in flatten(el):
                yield sub
        else:
            yield el

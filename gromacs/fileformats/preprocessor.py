# preprocessor.py, derived from pypreprocessor.py version 0.4.0
# http://code.google.com/p/pypreprocessor/
# Copyright (c) 2010 Evan Plaice
# Copyright (c) 2011 Oliver Beckstein
# Licence: MIT

"""
Preprocessor in Python
======================

:Author: Evan Plaice, Oliver Beckstein
:Licence: MIT
:Year: 2010, 2011

.. _pypreprocessor: http://code.google.com/p/pypreprocessor/

The :mod:`preprocessor` module is derived from `pypreprocessor`_ (release
0.4.0) and can be used to transform ITP file that make use of C-preprocessor
:ref:`directives-label` such as ::

  #ifdef POSRES
  [ position_restraints ]
  ...
  ...
  #endif

To include the position restraints use the :class:`Preprocessor`::

  PP = Preprocessor(filename="drug.itp", output="pp_drug.itp", POSRES=True)
  PP.parse()

To exclude them::

  PP = Preprocessor(filename="drug.itp", output="pp_drug.itp", POSRES=False)
  PP.parse()

(or simply omit ``POSRES`` from the argument list).

Once the file has been parsed, one can (1) write it out with the
:meth:`Preprocessor.write` method::

 PP.write("parsed.itp")

Or (2) create a :func:`cStringIO.StringIO` object from the in-memory parsed file
with :meth:`Preprocessor.StringIO`::::

 s = PP.StringIO()
 for line in s:
     print s,
 s.close()

Finally, there's also a `context manager`_ for the :mod:`cStringIO`
functionality, provided by :meth:`Preprocessor.open`::

 with PP.open() as s:
    for line in s:
        print line,

.. _`context manager`:
   http://docs.python.org/reference/datamodel.html#context-managers


.. _directives-label:

Directives
----------

The directives understood are:

``#define VAR``

  Define the variable *VAR*; note that even ``#define VAR 0`` in the
  input file will have the effect of defining the variable. The
  :class:`Preprocessor` constructor, however, will *not define* any
  *VAR* keyword that evaluates to ``False``.

``#undef VAR``

  Undefines *VAR*.

``#ifdef VAR`` ... ``#else`` ... ``#endif``

  Conditional evaluation of content blocks. *VAR* can only be a simple
  variable name, and it is checked against a list of defined variable
  names.

``#exclude`` ... ``#endexclude``

  Content inside the exclude block is omitted from the processed file.

.. Note::

   Expansion of a defined VAR inside the file is *not supported*. VAR are *only used
   in the context of ``#ifdef``/``#ifndef`` directives.


Classes
-------

.. autoclass:: Preprocessor
   :members:

"""
from __future__ import with_statement

__author__ = 'Evan Plaice'
__version__ = '0.4.0'
__licence__ = "MIT"
__URL__ = "http://code.google.com/p/pypreprocessor/"

import os
from contextlib import contextmanager

class Preprocessor(object):
    """CPP-style processing of files.

    The directives understood are:

    - ``#define VAR``
    - ``#undef VAR``
    - ``#ifdef VAR`` ... ``#else`` ... ``#endif``
    - ``#exclude`` ... ``#endexclude``

    """
    def __init__(self, filename, output=None, **kwargs):
        """Set up the preprocessor

        :Arguments:
           *filename*
              input file name
        :Keywords:
           *output*
              output file name, if ``None`` then a temporary file will be
              created for :meth:`write`.
           *clean*
              leave at ``True`` (to remove any preprocessor directives and
              anything excluded through #ifdefs); use ``False`` for debugging,
              which simply prefixes excluded lines with "*commentchar* #"
           *commentchar*
              how to comment out lines, leave at the default for itp files [";"]
           *strip*
              remove all empty lines and lines starting with *commentchar*
              (does not work with *clean* = ``False``) [``False``]
           *defines*
              any other keywords *VAR* are interpreted as ``#define VAR`` statement if
              *VAR* evaluates to ``True``.

        .. versionchanged:: 0.3.1
           *strip* keyword added
        """
        # public variables
        self.input = filename  # XXX: was a filename
        self.output = output
        self.removeMeta = kwargs.pop("clean", True)
        self.commentchar = kwargs.pop("commentchar", ";")  # for itp files
        self.default_defines = [x for x in kwargs if kwargs[x]]    #   #define x
        self.defines = self.default_defines[:]                     # see parser()
        self.strip = kwargs.pop('strip', False)
        if not self.removeMeta and self.strip:
            import warnings
            warnings.warn("Preprocessor: clean=False takes precedence over strip=True")
        # private variables
        self.__linenum = 0
        self.__excludeblock = False
        self.__ifblock = False
        self.__ifcondition = ''
        self.__ifconditions = []
        self.__evalsquelch = True
        self.__outputBuffer = ''

    def define(self, define):
        """#define directive"""
        self.defines.append(define)

    def search_defines(self, define):
        """Check if variable *define* has been defined."""
        return (define in self.defines)

    def compare_defines_and_conditions(self, defines, conditions):
        """#ifdef directive"""
        # if defines and conditions lists have no intersecting values (ie. else = true)
        return (not [val for val in defines if val in conditions])

    def undefine(self, define):
        """#undef directive"""
        # re-map the defines list excluding the define specified in the args
        self.defines[:] = [x for x in self.defines if x != define]

    def lexer(self, line):
        """evaluate *line*

        :returns: ``(squelch, metadata)``
        """
        if self.__ifblock is False and self.__excludeblock is False:
            if line[:1] != '#':
                return False, False
        # handle #define directives
        if line[:7] == '#define':
            if len(line.split()) != 2:
                self.raise_error('#define')
            else:
                self.define(line.split()[1])
                return False, True
        # handle #undef directives
        if line[:6] == '#undef':
            if len(line.split()) != 2:
                self.raise_error('#undef')
            else:
                self.undefine(line.split()[1])
                return False, True
        # handle #endif directives
        if line[:6] == '#endif':
            if len(line.split()) != 1:
                self.raise_error('#endif')
            else:
                self.__ifblock = False
                self.__ifcondition = ''
                self.__ifconditions = []
                return False, True
        # handle #endexclude directives
        if line[:11] == '#endexclude':
            if len(line.split()) != 1:
                self.raise_error('#endexclude')
            else:
                self.__excludeblock = False
                return False, True
        # handle #exclude directives
        if line[:8] == '#exclude':
            if len(line.split()) != 1:
                self.raise_error('#exclude')
            else:
                self.__excludeblock = True
        # process the excludeblock
        if self.__excludeblock is True:
            return True, False
        # handle #ifdef directives
        if line[:6] == '#ifdef':
            if len(line.split()) != 2:
                self.raise_error('#ifdef')
            else:
                self.__ifblock = True
                self.__ifcondition = line.split()[1]
                self.__ifconditions.append(line.split()[1])
        # handle #else directives
        if line[:5] == '#else':
            if len(line.split()) != 1:
                self.raise_error('#else')
        # process the ifblock
        if self.__ifblock is True:
            # evaluate and process an #ifdef
            if line[:6] == '#ifdef':
                self.__evalsquelch = not self.search_defines(self.__ifcondition)
                return False, True
            # evaluate and process the #else
            elif line[:5] == '#else':
                self.__evalsquelch = not self.compare_defines_and_conditions(self.defines, self.__ifconditions)
                return False, True
            else:
                return self.__evalsquelch, False
        else:
            return False, False

    def raise_error(self, directive):
        """error handling

        :Raises: :exc:`SyntaxError`
        """
        msg = 'File: "' + self.input + '", line ' + str(self.__linenum)
        msg = msg + '\n' + 'SyntaxError: Invalid ' + directive + ' directive'
        raise SyntaxError(msg)

    def parse(self, **kwargs):
        """parsing/processing

        *kwargs* are variables that are set (``#define VAR``) or unset
        (``#undef VAR``) at the beginning of the file. They are applied to all
        the defines that were provided to the constructor, and hence using
        *VAR* = ``False`` allows one to undefine some of these.

        This method only populates the output buffer and does not write an
        output file; use :meth:`write` for that purpose.
        """
        # unset defaults for any VAR=False
        self.defines = [x for x in self.default_defines if kwargs.pop(x, True)]
        # add all new defines for which VAR is True
        self.defines.extend([x for x in kwargs if kwargs[x]])

        self.__outputBuffer = ''
        with open(self.input, 'r') as input_file:
            # process the input file
            for line in input_file:
                self.__linenum += 1
                # to squelch or not to squelch
                squelch, metaData = self.lexer(line)
                # process and output
                if self.removeMeta is True:
                    if metaData is True or squelch is True:
                        continue
                if squelch is True:
                    self.__outputBuffer += self.commentchar + '#' + line
                    continue
                if self.strip and (len(line.strip()) == 0 or line.strip().startswith(self.commentchar)):
                    continue
                # output survived!
                self.__outputBuffer += line

    @property
    def buffer(self):
        """String representation of the processed input file.

        .. SeeAlso:: :meth:`StringIO` and :meth:`write`
        """
        return str(self.__outputBuffer)

    def StringIO(self):
        """Return a :func:`cStringIO.StringIO` instance of the buffer.

        The resulting instance can be treated as a read-only file containing
        the processed input from the last invocation of :meth:`parse`.
        """
        from cStringIO import StringIO
        return StringIO(self.__outputBuffer)

    @contextmanager
    def open(self):
        """Open the processed file (in memory) for reading.

        Use as ::

           with preprocessor.open() as pp:
              for line in pp:
                  print line,

        The method returns a :func:`cStringIO.StringIO` stream as provided by
        :meth:`StringIO`.
        """
        stream = self.StringIO()
        try:
            yield stream
        finally:
            stream.close()

    def write(self, filename=None):
        """write out file

        If *filename* is ``None`` then the constructor default (*output*) is
        chosen. If *filename* is ``False`` or no output filename was set then a
        temporary file is created (and it is the user's responsibility to clean
        up this file).

        Returns the filename of the file written.
        """
        if filename is None:
            filename = self.output

        if filename:
            output_file = open(filename, 'w')
        else:
            import tempfile
            root, ext = os.path.splitext(self.input)
            fd, filename = tempfile.mkstemp(suffix=ext, prefix="pp_", text=True)
            output_file = os.fdopen(fd, 'w')

        try:
            output_file.write(self.__outputBuffer)
        finally:
            output_file.close()

        return filename

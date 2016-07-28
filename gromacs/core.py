# GromacsWrapper: core.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
:mod:`gromacs.core` -- Core functionality
=========================================

Here the basic command class :class:`GromacsCommand` is defined. All Gromacs
command classes in :mod:`gromacs.tools` are automatically generated from
it. The documentation of :class:`GromacsCommand` applies to all wrapped Gromacs
commands and should be read by anyone using this package.


.. _input-output-label:

Input and Output
----------------

Each command wrapped by either :class:`GromacsCommand` or :class:`Command`
takes three additional keyword arguments: *stdout*, *stderr*, and
*input*. *stdout* and *stderr* determine how the command returns its own
output.

The *input* keyword is a string that is fed to the standard input of the
command (actually, :attr:`subprocess.Popen.stdin`). Or, if it is not string-like
then we assume it's actually a file-like object that we can read from, e.g. a
:attr:`subprocess.Popen.stdout` or a :class:`File`.

By setting the *stdout* and *stderr* keywords appropriately, one can have the
output simply printed to the screen (use ``True``; this is the default,
although see below for the use of the ``capture_output``
:mod:`gromacs.environment` flag), capture in a python variable as a string for
further processing (use ``False``), write to a file (use a :class:`File`
instance) or as input for another command (e.g. use the
:attr:`subprocess.Popen.stdin`).

When writing setup- and analysis pipelines it can be rather cumbersome to have
the gromacs output on the screen. For these cases GromacsWrapper allows you to
change its behaviour globally. By setting the value of the
:mod:`gromacs.environment` :class:`~gromacs.environment.Flag`
``capture_output`` to ``True`` (in the GromacsWrapper
:data:`gromacs.environment.flags` registry) ::

  import gromacs.environment
  gromacs.environment.flags['capture_output'] = True

all commands will capture their output (like *stderr* = ``False`` and *stdout*
= ``False``). Explicitly setting these keywords overrides the global
default. The default value for ``flags['capture_output']`` is ``False``,
i.e. output is directed through STDOUT and STDERR.

.. Warning::

   One downside of ``flags['capture_output'] = True`` is that it becomes much
   harder to debug scripts unless the script is written in such a way to show
   the output when the command fails. Therefore, it is advisable to only
   capture output on well-tested scripts.

A third value of ``capture_output`` is the value ``"file"``::

    gromacs.environment.flags['capture_output'] = "file"

This writes the captured output to a file. The file name is specified in
``flags['capture_output_filename'`` and defaults to
*"gromacs_captured_output.txt"*. This file is *over-written* for each
command. In this way one can investigate the output from the last command
(presumably because it failed). STDOUT and STDERR are captured into this file
by default. STDERR is printed first and then STDOUT, which does not necessarily
reflect the order of output one would see on the screen. If your code captures
STDOUT for further processing then an uncaptured STDERR is written to the
capture file.

.. Note::

   There are some commands for which capturing output
   (``flags['capture_output'] = True``) might be problematic. If the command
   produces a large or inifinite amount of data then a memory error will occur
   because Python nevertheless stores the output internally first. Thus one
   should avoid capturing progress output from
   e.g. :class:`~gromacs.tools.Mdrun` unless the output has been throttled
   appropriately.


Classes
-------

.. autoclass:: GromacsCommand
   :members: __call__, run, transform_args, Popen, help,
             check_failure, gmxdoc
   :inherited-members:

.. autoclass:: Command
   :members:  __call__, run, transform_args, Popen, help,
             command_name

.. autoclass:: PopenWithInput
   :members:
"""
from __future__ import absolute_import, with_statement

__docformat__ = "restructuredtext en"


import itertools
import logging
import re
import subprocess

from subprocess import PIPE
from .exceptions import GromacsError


DEFAULT_OUTPUT = None
DEFAULT_ERROR_OUTPUT = None

RAISE = 1
WARN = 2
IGNORE = 3

FAILURE_MODE = RAISE


def set_default_output(out):
    global DEFAULT_OUTPUT
    DEFAULT_OUTPUT = out


def set_default_error_output(out):
    global DEFAULT_ERROR_OUTPUT
    DEFAULT_ERROR_OUTPUT = out


def set_failure_mode(mode):
    global FAILURE_MODE
    FAILURE_MODE = mode


logger = logging.getLogger(__name__)


class Command(object):

    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs
        self.last_command = None
        self.last_output = None
        self.last_error = None
        self.last_returncode = None

    def run(self, *args, **kwargs):
        stdin, input = self._get_stdin(kwargs.pop('input', None))
        stdout = kwargs.pop('stdout', DEFAULT_OUTPUT)
        stderr = kwargs.pop('stderr', DEFAULT_ERROR_OUTPUT)
        cmd = self.last_command = self.make_command_line(args, kwargs)

        proc = subprocess.Popen(cmd, stdin=stdin, stdout=stdout, stderr=stderr)
        output, err = proc.communicate(input=input)

        self.last_output = output
        self.last_error = err
        self.last_returncode = proc.returncode

        return proc.returncode, output, err

    def __call__(self, *args, **kwargs):
        self.run(*args, **kwargs)

    def make_command_line(self, args, kwargs):
        cmd = []

        for arg in self.args + args:
            cmd.extend(self.expand_param(arg))
        for opt, value in self.kwargs.items():
            cmd.extend(['-%s' % opt] + self.expand_param(value))
        for opt, value in kwargs.items():
            cmd.extend(['-%s' % opt] + self.expand_param(value))

        return cmd

    def check_failure(self):
        if self.last_returncode != 0:
            output = [self.last_output, self.last_error]
            output = '\n'.join([x for x in output if x is not None])

            last_cmd = ' '.join(self.last_command)
            message = '%s --> %d' % (last_cmd, self.last_returncode)

            if FAILURE_MODE == IGNORE:
                pass
            if FAILURE_MODE == WARN:
                logger.warn(message)
            else:
                raise GromacsError(message)


    @classmethod
    def expand_param(cls, param):
        if type(param) is bool:
            return [str(param).lower()]
        elif type(param) in [list, tuple]:
            params = [cls.expand_param(p) for p in param]
            return itertools.chain.from_iterable(params)
        else:
            return [str(param)]

    @classmethod
    def _get_stdin(cls, input):
        if input:
            if isinstance(input, basestring):
                return PIPE, input
            try:
                return PIPE, '\n'.join(input)
            except TypeError:
                raise ValueError("input must be a string, None or an iterable")
        else:
            return None, None


class GromacsCommand(Command):

    command_name = None
    driver = None

    documentation_re = r".*?(?P<DOCS>DESCRIPTION.*)"
    failure_re = r"""----+\n                        # ---- decorator line
            \s*Program\s+(?P<program_name>\w+),     #  Program name,
              \s+VERSION\s+(?P<version>[\w.]+)\s*\n #    VERSION 4.0.5
            (?P<message>.*?)\n                      # full message, multiple lines
            \s*                                     # empty line (?)
            ----+\n                                 # ---- decorator line
            """

    def __init__(self, *args, **kwargs):
        if self.driver:
            args = [self.driver, self.command_name] + list(args)
        else:
            args = [self.command_name] + list(args)
        super(GromacsCommand, self).__init__(*args, **kwargs)
        self._doc_cache = None

    def check_failure(self):
        if self.last_returncode != 0:
            output = [self.last_output, self.last_error]
            output = '\n'.join([x for x in output if x is not None])

            last_cmd = ' '.join(self.last_command)
            message = '%s --> %d' % (last_cmd, self.last_returncode)

            m = re.search(self.failure_re, output, re.VERBOSE | re.DOTALL)
            if m:
                message = '%s\n%s' % (message, m.group('message'))

            if FAILURE_MODE == IGNORE:
                pass
            if FAILURE_MODE == WARN:
                logger.warn(message)
            else:
                raise GromacsError(message)
        return self.last_returncode != 0

    def run(self, *args, **kwargs):
        result = super(GromacsCommand, self).run(*args, **kwargs)
        self.check_failure()
        return result

    def _get_documentation(self):
        if self._doc_cache is not None:
            return self._doc_cache
        rc, out, err = self.run(h=True, stdout=PIPE, stderr=PIPE)
        m = re.match(self.documentation_re, out, re.DOTALL)
        if not m:
            m = re.match(self.documentation_re, err, re.DOTALL)
        self._doc_cache = m.group('DOCS')
        return self._doc_cache
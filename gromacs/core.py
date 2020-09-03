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
from __future__ import absolute_import, with_statement, print_function
import six

__docformat__ = "restructuredtext en"

import sys
import re
import subprocess
from subprocess import STDOUT, PIPE
import warnings
import errno

import logging
logger = logging.getLogger('gromacs.core')


from .exceptions import GromacsError, GromacsFailureWarning
from . import environment

class Command(object):
    """Wrap simple script or command."""
    #: Derive a class from command; typically one only has to set *command_name*
    #: to the name of the script or executable. The full path is required if it
    #: cannot be found by searching :envvar:`PATH`.
    command_name = None

    def __init__(self, *args, **kwargs):
        """Set up the command class.

        The arguments can always be provided as standard positional
        arguments such as

          ``"-c", "config.conf", "-o", "output.dat", "--repeats=3", "-v", "input.dat"``

        In addition one can also use keyword arguments such as

          ``c="config.conf", o="output.dat", repeats=3, v=True``

        These are automatically transformed appropriately according to
        simple rules:

        * Any single-character keywords are assumed to be POSIX-style
          options and will be prefixed with a single dash and the value
          separated by a space.

        * Any other keyword is assumed to be a GNU-style long option
          and thus will be prefixed with two dashes and the value will
          be joined directly with an equals sign and no space.

        If this does not work (as for instance for the options of the
        UNIX ``find`` command) then provide options and values in the
        sequence of positional arguments.


        *Example*

        Create a ``Ls`` class whose instances execute the ``ls`` command::

          LS = type("LS", (gromacs.core.Command,), {'command_name': 'ls'})
          ls = LS()
          ls()        # lists directory like ls
          ls(l=True)  # lists directory like ls -l

        Now create an instance that performs a long directory listing by
        default::

          lslong = LS(l=True)
          lslong()    # like ls -l

        """

        self.args = args
        self.kwargs = kwargs

    def run(self, *args, **kwargs):
        """Run the command; args/kwargs are added or replace the ones given to the constructor."""
        _args, _kwargs = self._combine_arglist(args, kwargs)
        results, p = self._run_command(*_args, **_kwargs)
        return results

    def _combine_arglist(self, args, kwargs):
        """Combine the default values and the supplied values."""
        _args = self.args + args
        _kwargs = self.kwargs.copy()
        _kwargs.update(kwargs)
        return _args, _kwargs

    def _run_command(self, *args, **kwargs):
        """Execute the command; see the docs for __call__.

        :Returns: a tuple of the *results* tuple ``(rc, stdout, stderr)`` and
                  the :class:`Popen` instance.
        """
        # hack to run command WITHOUT input (-h...) even though user defined
        # input (should have named it "ignore_input" with opposite values...)
        use_input = kwargs.pop('use_input', True)

        # logic for capturing output (see docs on I/O and the flags)
        capturefile = None
        if environment.flags['capture_output'] is True:
            # capture into Python vars (see subprocess.Popen.communicate())
            kwargs.setdefault('stderr', PIPE)
            kwargs.setdefault('stdout', PIPE)
        elif environment.flags['capture_output'] == "file":
            if 'stdout' in kwargs and 'stderr' in kwargs:
                pass
            else:
                # XXX: not race or thread proof; potentially many commands write to the same file
                fn = environment.flags['capture_output_filename']
                capturefile = open(fn, "w")   # overwrite (clobber) capture file
                if 'stdout' in kwargs and 'stderr' not in kwargs:
                    # special case of stdout used by code but stderr should be captured to file
                    kwargs.setdefault('stderr', capturefile)
                else:
                    # merge stderr with stdout and write stdout to file
                    # (stderr comes *before* stdout in capture file, could split...)
                    kwargs.setdefault('stderr', STDOUT)
                    kwargs.setdefault('stdout', capturefile)

        try:
            p = self.Popen(*args, **kwargs)
            out, err = p.communicate(use_input=use_input) # special Popen knows input!
        except:
            if capturefile is not None:
                logger.error("Use captured command output in %r for diagnosis.", capturefile)
            raise
        finally:
            if capturefile is not None:
                capturefile.close()
        rc = p.returncode
        return (rc, out, err), p

    def _commandline(self, *args, **kwargs):
        """Returns the command line (without pipes) as a list."""
         # transform_args() is a hook (used in GromacsCommand very differently!)
        return [self.command_name] + self.transform_args(*args, **kwargs)

    def commandline(self, *args, **kwargs):
        """Returns the commandline that run() uses (without pipes)."""
        # this mirrors the setup in run()
        _args, _kwargs = self._combine_arglist(args, kwargs)
        return self._commandline(*_args, **_kwargs)

    def Popen(self, *args, **kwargs):
        """Returns a special Popen instance (:class:`PopenWithInput`).

        The instance has its input pre-set so that calls to
        :meth:`~PopenWithInput.communicate` will not need to supply
        input. This is necessary if one wants to chain the output from
        one command to an input from another.

        :TODO:
          Write example.
        """
        stderr = kwargs.pop('stderr', None)     # default: print to stderr (if STDOUT then merge)
        if stderr is False:                     # False: capture it
            stderr = PIPE
        elif stderr is True:
            stderr = None                       # use stderr

        stdout = kwargs.pop('stdout', None)     # either set to PIPE for capturing output
        if stdout is False:                     # ... or to False
            stdout = PIPE
        elif stdout is True:
            stdout = None                       # for consistency, make True write to screen

        stdin = kwargs.pop('stdin', None)
        input = kwargs.pop('input', None)

        use_shell = kwargs.pop('use_shell', False)
        if input:
            stdin = PIPE
            if isinstance(input, six.string_types) and not input.endswith('\n'):
                # make sure that input is a simple string with \n line endings
                input = six.text_type(input) + '\n'
            else:
                try:
                    # make sure that input is a simple string with \n line endings
                    input = '\n'.join(map(six.text_type, input)) + '\n'
                except TypeError:
                    # so maybe we are a file or something ... and hope for the best
                    pass

        cmd = self._commandline(*args, **kwargs)   # lots of magic happening here
                                                   # (cannot move out of method because filtering of stdin etc)
        try:
            p = PopenWithInput(cmd, stdin=stdin, stderr=stderr, stdout=stdout,
                               universal_newlines=True, input=input, shell=use_shell)
        except OSError as err:
            logger.error(" ".join(cmd))            # log command line
            if err.errno == errno.ENOENT:
                errmsg = "Failed to find Gromacs command {0!r}, maybe its not on PATH or GMXRC must be sourced?".format(self.command_name)
                logger.fatal(errmsg)
                raise OSError(errmsg)
            else:
                logger.exception("Setting up Gromacs command {0!r} raised an exception.".format(self.command_name))
                raise
        logger.debug(p.command_string)
        return p

    def transform_args(self, *args, **kwargs):
        """Transform arguments and return them as a list suitable for Popen."""
        options = []
        for option,value in kwargs.items():
            if not option.startswith('-'):
                # heuristic for turning key=val pairs into options
                # (fails for commands such as 'find' -- then just use args)
                if len(option) == 1:
                    option = '-' + option         # POSIX style
                else:
                    option = '--' + option        # GNU option
            if value is True:
                options.append(option)
                continue
            elif value is False:
                raise ValueError('A False value is ambiguous for option {0!r}'.format(option))

            if option[:2] == '--':
                options.append(option + '=' + str(value))    # GNU option
            else:
                options.extend((option, str(value)))         # POSIX style
        return options + list(args)

    def help(self, long=False):
        """Print help; same as using ``?`` in ``ipython``. long=True also gives call signature."""
        print("\ncommand: {0!s}\n\n".format(self.command_name))
        print(self.__doc__)
        if long:
            print("\ncall method: command():\n")
            print(self.__call__.__doc__)

    def __call__(self,*args,**kwargs):
        """Run command with the given arguments::

           rc,stdout,stderr = command(*args, input=None, **kwargs)

        All positional parameters *args* and all gromacs *kwargs* are passed on
        to the Gromacs command. input and output keywords allow communication
        with the process via the python subprocess module.

        :Arguments:
          *input* : string, sequence
             to be fed to the process' standard input;
             elements of a sequence are concatenated with
             newlines, including a trailing one    [``None``]
          *stdin*
             ``None`` or automatically set to ``PIPE`` if input given [``None``]
          *stdout*
             how to handle the program's stdout stream [``None``]

             filehandle
                    anything that behaves like a file object
             ``None`` or ``True``
                    to see  output on screen
             ``False`` or ``PIPE``
                     returns the output as a string in  the stdout parameter

          *stderr*
             how to handle the stderr stream [``None``]

             ``STDOUT``
                     merges standard error with the standard out stream
             ``False`` or ``PIPE``
                     returns the output as a string in the stderr return parameter
             ``None`` or ``True``
                     keeps it on stderr (and presumably on screen)

        Depending on the value of the GromacsWrapper flag
        :data:`gromacs.environment.flags```['capture_output']`` the above
        default behaviour can be different.

        All other kwargs are passed on to the Gromacs tool.

        :Returns:

           The shell return code rc of the command is always returned. Depending
           on the value of output, various strings are filled with output from the
           command.

        :Notes:

           In order to chain different commands via pipes one must use the special
           :class:`PopenWithInput` object (see :meth:`GromacsCommand.Popen` method) instead of the simple
           call described here and first construct the pipeline explicitly and then
           call the :meth:`PopenWithInput.communicate` method.

           ``STDOUT`` and ``PIPE`` are objects provided by the :mod:`subprocess` module. Any
           python stream can be provided and manipulated. This allows for chaining
           of commands. Use ::

              from subprocess import PIPE, STDOUT

           when requiring these special streams (and the special boolean
           switches ``True``/``False`` cannot do what you need.)

           (TODO: example for chaining commands)
        """
        return self.run(*args, **kwargs)


class GromacsCommand(Command):
    """Base class for wrapping a Gromacs tool.

    Limitations: User must have sourced ``GMXRC`` so that the python script can
    inherit the environment and find the gromacs programs.

    The class doc string is dynamically replaced by the documentation of the
    gromacs command the first time the doc string is requested. If the tool is
    not available at the time (i.e., cannot be found on :env:`PATH`) then the
    generic doc string is shown and an :exc:`OSError` exception is only raised
    when the user is actually trying to the execute the command.
    """

    # TODO: setup the environment from GMXRC (can use env=DICT in Popen/call)

    command_name = None
    driver = None
    doc_pattern = """.*?(?P<DOCS>DESCRIPTION.*)"""
    gmxfatal_pattern = """----+\n                   # ---- decorator line
            \s*Program\s+(?P<program_name>\w+),     #  Program name,
              \s+VERSION\s+(?P<version>[\w.]+)\s*\n #    VERSION 4.0.5
            (?P<message>.*?)\n                      # full message, multiple lines
            \s*                                     # empty line (?)
            ----+\n                                 # ---- decorator line
            """
    # matches gmx_fatal() output
    # -------------------------------------------------------
    # Program <program_name>, VERSION <version>
    # ... <message>
    # -------------------------------------------------------

    #: Available failure modes.
    failuremodes = ('raise', 'warn', None)

    def __init__(self, *args, **kwargs):
        """Set up the command with gromacs flags as keyword arguments.

        The following  are generic instructions; refer  to the Gromacs
        command  usage information  that should  have  appeared before
        this generic documentation.

        As an example, a generic Gromacs command could use the following flags::

          cmd = GromacsCommand('v', f=['md1.xtc','md2.xtc'], o='processed.xtc', t=200, ...)

        which would correspond to running the command in the shell as ::

          GromacsCommand -v -f md1.xtc md2.xtc -o processed.xtc -t 200

        **Gromacs command line arguments**

           Gromacs boolean switches (such as ``-v``) are given as python
           positional arguments (``'v'``) or as keyword argument (``v=True``);
           note the quotes in the first case. Negating a boolean switch can be
           done with ``'nov'``, ``nov=True`` or ``v=False`` (and even ``nov=False``
           works as expected: it is the same as ``v=True``).

           Any Gromacs options that take parameters are handled as keyword
           arguments. If an option takes multiple arguments (such as the
           multi-file input ``-f file1 file2 ...``) then the list of files must be
           supplied as a python list.

           If a keyword has the python value ``None`` then it will *not* be
           added to the Gromacs command line; this allows for flexible
           scripting if it is not known in advance if an input file is
           needed. In this case the default value of the gromacs tool
           is used.

           Keywords must be legal python keywords or the interpreter raises a
           :exc:`SyntaxError` but of course Gromacs commandline arguments are
           not required to be legal python. In this case "quote" the option
           with an underscore (``_``) and the underscore will be silently
           stripped. For instance, ``-or`` translates to the illegal keyword
           ``or`` so it must be underscore-quoted::

              cmd(...., _or='mindistres.xvg')

        **Command execution**

           The command is executed with the :meth:`~GromacsCommand.run` method or by
           calling it as a function. The two next lines are equivalent::

             cmd(...)
             cmd.run(...)

           When the command is run one can override options that were given at
           initialization or one can add additional ones. The same rules for
           supplying Gromacs flags apply as described above.

        **Non-Gromacs keyword arguments**

           The other keyword arguments (listed below) are not passed on to the
           Gromacs tool but determine how the command class behaves. *They are
           only useful when instantiating a class*, i.e. they determine how
           this tool behaves during all future invocations although it can be
           changed by setting :attr:`failuremode`. This is mostly of interest
           to developers.

        :Keywords:
           *failure*
              determines how a failure of the gromacs command is treated; it
              can be one of the following:

              'raise'
                   raises GromacsError if command fails
              'warn'
                   issue a :exc:`GromacsFailureWarning`
              ``None``
                   just continue silently

           *doc* : string
              additional documentation (*ignored*) []

        .. versionchanged:: 0.6.0
           The *doc* keyword is now ignored (because it was not worth the effort to
           make it work with the lazy-loading of docs).
        """
        doc = kwargs.pop('doc', None)  # ignored
        self.__failuremode = None
        self.failuremode = kwargs.pop('failure', 'raise')
        self.gmxargs = self._combineargs(*args, **kwargs)
        self._doc_cache = None

    def failuremode():
        doc = """mode determines how the GromacsCommand behaves during failure

        It can be one of the following:

              'raise'
                   raises GromacsError if command fails
              'warn'
                   issue a :exc:`GromacsFailureWarning`
              ``None``
                   just continue silently

        """
        def fget(self):
            return self.__failuremode
        def fset(self, mode):
            if not mode in self.failuremodes:
                raise ValueError('failuremode must be one of {0!r}'.format(self.failuremodes))
            self.__failuremode = mode
        return locals()
    failuremode = property(**failuremode())

    def _combine_arglist(self, args, kwargs):
        """Combine the default values and the supplied values."""
        gmxargs = self.gmxargs.copy()
        gmxargs.update(self._combineargs(*args, **kwargs))
        return (), gmxargs    # Gromacs tools don't have positional args --> args = ()

    def check_failure(self, result, msg='Gromacs tool failed', command_string=None):
        rc, out, err = result
        if command_string is not None:
            msg += '\nCommand invocation: ' + str(command_string)
        had_success = (rc == 0)
        if not had_success:
            gmxoutput = "\n".join([x for x in [out, err] if x is not None])
            m = re.search(self.gmxfatal_pattern, gmxoutput, re.VERBOSE | re.DOTALL)
            if m:
                formatted_message = ['GMX_FATAL  '+line for line in m.group('message').split('\n')]
                msg = "\n".join(\
                    [msg, "Gromacs command {program_name!r} fatal error message:".format(**m.groupdict())] +
                    formatted_message)
            if self.failuremode == 'raise':
                raise GromacsError(rc, msg)
            elif self.failuremode == 'warn':
                warnings.warn(msg + '\nError code: {0!r}\n'.format(rc), category=GromacsFailureWarning)
            elif self.failuremode is None:
                pass
            else:
                raise ValueError('unknown failure mode {0!r}'.format(self.failuremode))
        return had_success

    def _combineargs(self, *args, **kwargs):
        """Add switches as 'options' with value True to the options dict."""
        d = {arg: True for arg in args}   # switches are kwargs with value True
        d.update(kwargs)
        return d

    def _build_arg_list(self, **kwargs):
        """Build list of arguments from the dict; keys must be valid  gromacs flags."""
        arglist = []
        for flag, value in kwargs.items():
            # XXX: check flag against allowed values
            flag = str(flag)
            if flag.startswith('_'):
                flag = flag[1:]                 # python-illegal keywords are '_'-quoted
            if not flag.startswith('-'):
                flag = '-' + flag               # now flag is guaranteed to start with '-'
            if value is True:
                arglist.append(flag)            # simple command line flag
            elif value is False:
                if flag.startswith('-no'):
                    # negate a negated flag ('noX=False' --> X=True --> -X ... but who uses that?)
                    arglist.append('-' + flag[3:])
                else:
                    arglist.append('-no' + flag[1:])  # gromacs switches booleans by prefixing 'no'
            elif value is None:
                pass                            # ignore flag = None
            else:
                try:
                    arglist.extend([flag] + value) # option with value list
                except TypeError:
                    arglist.extend([flag, value])  # option with single value
        return list(map(str, arglist))  # all arguments MUST be strings

    def _run_command(self,*args,**kwargs):
        """Execute the gromacs command; see the docs for __call__."""
        result, p = super(GromacsCommand, self)._run_command(*args, **kwargs)
        self.check_failure(result, command_string=p.command_string)
        return result, p

    def _commandline(self, *args, **kwargs):
        """Returns the command line (without pipes) as a list. Inserts driver if present"""
        if(self.driver is not None):
            return [self.driver, self.command_name] + self.transform_args(*args, **kwargs)
        return [self.command_name] + self.transform_args(*args, **kwargs)


    def transform_args(self,*args,**kwargs):
        """Combine arguments and turn them into gromacs tool arguments."""
        newargs = self._combineargs(*args, **kwargs)
        return self._build_arg_list(**newargs)

    def _get_gmx_docs(self):
        """Extract standard gromacs doc

        Extract by running the program and chopping the header to keep from
        'DESCRIPTION' onwards.
        """
        if self._doc_cache is not None:
            return self._doc_cache

        try:
            logging.disable(logging.CRITICAL)
            rc, header, docs = self.run('h', stdout=PIPE, stderr=PIPE, use_input=False)
        except:
            logging.critical("Invoking command {0} failed when determining its doc string. Proceed with caution".format(self.command_name))
            self._doc_cache = "(No Gromacs documentation available)"
            return self._doc_cache
        finally:
            # ALWAYS restore logging...
            logging.disable(logging.NOTSET)

        # The header is on STDOUT and is ignored. The docs are read from STDERR in GMX 4.
        m = re.match(self.doc_pattern, docs, re.DOTALL)

        if m is None:
            # In GMX 5, the opposite is true (Grrr)
            m = re.match(self.doc_pattern, header, re.DOTALL)
            if m is None:
                self._doc_cache = "(No Gromacs documentation available)"
                return self._doc_cache

        self._doc_cache = m.group('DOCS')
        return self._doc_cache


class PopenWithInput(subprocess.Popen):
    """Popen class that knows its input.

    1. Set up the instance, including all the input it shoould receive.
    2. Call :meth:`PopenWithInput.communicate` later.

    .. Note:: Some versions of python have a bug in the subprocess module
              ( `issue 5179`_ ) which does not clean up open file
              descriptors. Eventually code (such as this one) fails with the
              error:

                  *OSError: [Errno 24] Too many open files*

              A weak workaround is to increase the available number of open
              file descriptors with ``ulimit -n 2048`` and run analysis in
              different scripts.

    .. _issue 5179: http://bugs.python.org/issue5179
    """

    def __init__(self, *args, **kwargs):
        """Initialize with the standard :class:`subprocess.Popen` arguments.

        :Keywords:
           *input*
               string that is piped into the command

        """
        kwargs.setdefault('close_fds', True)   # fixes 'Too many open fds' with 2.6
        self.input = kwargs.pop('input', None)
        if six.PY2 and self.input is not None:
            # in Python 2, subprocess.Popen uses os.write(chunk) with default ASCII encoding
            self.input = self.input.encode('utf-8')
        self.command = args[0]
        try:
            input_string = 'printf "' + \
                self.input.replace('\n','\\n') + '" | '  # display newlines
        except (TypeError, AttributeError):
            input_string = ""
        self.command_string = input_string + " ".join(self.command)
        super(PopenWithInput,self).__init__(*args, **kwargs)

    def communicate(self, use_input=True):
        """Run the command, using the input that was set up on __init__ (for *use_input* = ``True``)"""
        if use_input:
            return super(PopenWithInput, self).communicate(self.input)
        else:
            return super(PopenWithInput, self).communicate()

    def __str__(self):
        return "<Popen on {0!r}>".format(self.command_string)

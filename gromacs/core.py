# $Id$
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
:mod:`gromacs.core` -- Core functionality
======================================

Here the basic command class :class:`GromacsCommand` is defined. All
Gromacs command classes in :mod:`gromacs.tools` are automatically
generated from it.

.. autoclass:: GromacsCommand

"""
__docformat__ = "restructuredtext en"

import sys
import re
import subprocess
from subprocess import STDOUT, PIPE
import warnings
import errno

from gromacs import GromacsError, GromacsFailureWarning

class GromacsCommand(object):
    """Base class for wrapping a g_* command.
    
    Limitations: User must have sourced ``GMXRC`` so that the python script can
    inherit the environment and find the gromacs programs.

    The class doc string is dynamically replaced by the documentation of the
    gromacs command when an instance is created.
    """
    # TODO: setup the environment from GMXRC (can use env=DICT in Popen/call)

    command_name = None
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

    failuremodes = ('raise', 'warn', None)

    def __init__(self,*args,**kwargs):
        """Set up the command with gromacs flags as keyword arguments::

          cmd = GromacsCommand('v', f=['md1.xtc','md2.xtc'], o='processed.xtc', t=200, ...)

        Gromacs command line arguments

           Gromacs boolean switches (such as ``-v``) are given as python
           positional arguments (``'v'``) or as keyword argument (``v=True``);
           note the quotes in the first case. Negating at boolean switch can be
           done with ``'-nov'``, ``nov=True`` or ``v=False``.

           Any Gromacs options that take parameters are handled as keyword
           arguments. If an option takes multiple arguments (such as the
           multi-file input ``-f file1 file2 ...``) then the list of files must be
           supplied as a python list.

           If a keyword has the python value None then it will *not* be added to
           the Gromacs command line; this allows for flexible scripting if it is
           not known in advance if an input file is needed.

        Command execution

           The command is executed with the run() method or by
           calling it as a function. The two next lines are equivalent::

             cmd(...)
             cmd.run(...)

           When the command is run one can override options that were given at
           initialization or add additional ones.

        Non-Gromacs keyword arguments

           The other keyword arguments are not passed on to the Gromacs tool
           but determine how the command class behaves.

                
        :Keywords:
           failure
              determines how a failure of the gromacs command is treated; it
              can be one of the following:

              'raise'
                   raises GromacsError if command fails
              'warn'
                   issue a ``GromacsFailureWarning``
              ``None``
                   just continue silently

           doc : string
              additional documentation []
        """

        self.failuremode = kwargs.pop('failure','raise')
        self.extra_doc = kwargs.pop('doc',None)
        if not self.failuremode in self.failuremodes:
            raise ValueError('failuremode must be one of\n%(failuremodes)r' % vars(self))
        self.gmxargs = self._combineargs(*args, **kwargs)
        self.__doc__ = self.gmxdoc

    def run(self,*args,**kwargs):
        """Run the command; kwargs are added or replace the ones given to the constructor."""
        gmxargs = self.gmxargs.copy()
        gmxargs.update(self._combineargs(*args,**kwargs))
        return self._run_command(**gmxargs)

    def check_failure(self, result, msg='Gromacs tool failed', command_string=None):
        rc, out, err = result
        if not command_string is None:
            msg += '\nCommand invocation: ' + str(command_string)
        had_success = (rc == 0)
        if not had_success:
            gmxoutput = "\n".join([x for x in [out, err] if not x is None])
            m = re.search(self.gmxfatal_pattern, gmxoutput, re.VERBOSE | re.DOTALL)
            if m:
                formatted_message = ['GMX_FATAL  '+line for line in m.group('message').split('\n')]
                msg = "\n".join(\
                    [msg, "Gromacs command %(program_name)r fatal error message:" % m.groupdict()] +
                    formatted_message)
            if self.failuremode == 'raise':                
                raise GromacsError(rc, msg)
            elif self.failuremode == 'warn':
                warnings.warn(msg + '\nError code: %r\n' % rc, category=GromacsFailureWarning)
            elif self.failuremode is None:
                pass
            else:
                raise ValueError('unknown failure mode %r' % self.failuremode)
        return had_success
            

    def _combineargs(self,*args,**kwargs):
        """Add switches as 'options' with value True to the options dict."""
        d = dict([(arg, True) for arg in args])   # switches are kwargs with value True
        d.update(kwargs)
        return d
    
    def _build_arg_list(self,**kwargs):
        """Build list of arguments from the dict; keys must be valid  gromacs flags."""
        arglist = []
        for flag,value in kwargs.items():
            # XXX: check flag against allowed values
            flag = str(flag)
            if not flag.startswith('-'):
                flag = '-' + flag
            if value is True:
                arglist.append(flag)            # simple command line flag
            elif value is False:
                # XXX: does not work for '-noXXX False' ... but who uses that?
                arglist.append('-no'+flag[1:])  # gromacs switches booleans by prefixing 'no'
            elif value is None:
                pass                            # ignore flag = None
            else:
                try:
                    arglist.extend([flag] + value) # option with value list
                except TypeError:
                    arglist.extend([flag, value])  # option with single value
        return map(str, arglist)  # all arguments MUST be strings 

    def _run_command(self,*args,**kwargs):
        """Execute the gromacs command; see the docs for __call__."""
        use_input = kwargs.pop('use_input', True)     # hack to run command WITHOUT input (-h...)
        p = self.Popen(*args, **kwargs)
        out, err = p.communicate(use_input=use_input) # special Popen knows input!
        rc = p.returncode
        result = rc, out, err
        self.check_failure(result, command_string=p.command_string)
        return result

    def Popen(self, *args,**kwargs):
        """Returns a special Popen instance with pre-set input for communicate()."""

        stderr = kwargs.pop('stderr', STDOUT)   # default: Merge with stdout
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
        if input:
            stdin = PIPE
            if type(input) is str:
                # make sure that input is a simple string with \n line endings
                if not input.endswith('\n'):
                    input += '\n'
            else:
                try:
                    # make sure that input is a simple string with \n line endings
                    input = '\n'.join(map(str, input)) + '\n'
                except TypeError:
                    # so maybe we are a file or something ... and hope for the best
                    pass

        newargs = self._combineargs(*args,**kwargs)
        cmd = [self.command_name] + self._build_arg_list(**newargs)
        try:
            p = PopenWithInput(cmd, stdin=stdin, stderr=stderr, stdout=stdout,
                               universal_newlines=True, input=input)
        except OSError,err:
            if err.errno == errno.ENOENT:
                raise OSError("Failed to find Gromacs command '%r'. Source GMXRC." %
                              self.command_name)
            else:
                raise
        return p
        

    def _get_gmx_docs(self):
        """Extract standard gromacs doc by running the program and chopping the header."""        
        # Uses the class-wide arguments so that 'canned invocations' in cbook
        # are accurately reflected. Might be a problem when these invocations
        # supply wrong arguments... TODO: maybe check rc for that?
        # use_input=False needed for running commands in cbook that have input pre-defined
        rc,docs,nothing = self.run('h', stdout=PIPE, use_input=False)
        m = re.match(self.doc_pattern, docs, re.DOTALL)    # keep from DESCRIPTION onwards
        if m is None:
            return "(No Gromacs documentation available)"
        return m.group('DOCS')
        

    def gmxdoc():
        doc = """Usage for the underlying Gromacs tool (cached)."""
        def fget(self):
            if not (hasattr(self, '__doc_cache') and self.__doc_cache):
                self.__doc_cache = self._get_gmx_docs()
            docs = self.__doc_cache
            if self.extra_doc:
                docs = '\n'.join([self.extra_doc,'',
                                  "Documentation of the gromacs tool", 34*'=',
                                  docs])
            return docs
        return locals()
    gmxdoc = property(**gmxdoc())

    def help(self,long=False):
        """Print help; same as using ``?`` in ``ipython``. long=True also gives call signature."""
        print "\ncommand: %s\n\n" % self.command_name
        print self.__doc__
        if long:
            print "\ncall method: command():\n"
            print self.__call__.__doc__
        
    def __call__(self,*args,**kwargs):
        """Run command with the given arguments::

           rc,stdout,stderr = command(*args, input=None, **kwargs)
           
        All positional parameters \*args and all gromacs \*\*kwargs are passed on
        to the Gromacs command. input and output keywords allow communication
        with the process via the python subprocess module.
        
        :Arguments:
          input : string, sequence            
             to be fed to the process' standard input;
             elements of a sequence are concatenated with
             newlines, including a trailing one    [``None``]
          stdin
             ``None`` or automatically set to ``PIPE`` if input given [``None``]
          stdout
             how to handle the program's stdout stream [``None``]

             filehandle
                    anything that behaves like a file object
             ``None`` or ``True``
                    to see  output on screen
             ``False`` or ``PIPE``
                     returns the output as a string in  the stdout parameter 

          stderr
             how to handle the stderr stream [``STDOUT``]

             ``STDOUT``
                     merges standard error with the standard out stream
             ``False`` or ``PIPE``
                     returns the output as a string in the stderr return parameter
             ``None`` or ``True``
                     keeps it on stderr (and presumably on screen)

        All other kwargs are passed on to the Gromacs tool.
     
        :Returns:

           The shell return code rc of the command is always returned. Depending
           on the value of output, various strings are filled with output from the
           command.

        :Notes:

           By default, the process stdout and stderr are merged.

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
        return self.run(*args,**kwargs)


class PopenWithInput(subprocess.Popen):
    """Popen class that knows its input; simply call communicate() later."""

    def __init__(self,*args,**kwargs):
        """Initialize with the standard :class:`subprocess.Popen` arguments and *input*."""
        self.input = kwargs.pop('input',None)
        self.command = args[0]
        self.command_string = " ".join(self.command)
        super(PopenWithInput,self).__init__(*args,**kwargs)
    def communicate(self, use_input=True):
        if use_input:
            return super(PopenWithInput,self).communicate(self.input)
        else:
            return super(PopenWithInput,self).communicate()
    def __str__(self):
        return "<Popen on %r>" % self.command_string

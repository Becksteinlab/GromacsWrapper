# $Id$
"""Core functionality for the Gromacs python shell."""

import re
import subprocess
from subprocess import STDOUT, PIPE
import warnings
import errno

class GromacsCommand(object):
    """Base class for wrapping a g_* command.
    
    Limitations: User must have sourced GMXRC so that the python script can
    inherit the environment and find the gromacs programs.

    The class doc string is dynamically replaced by the documentation of the
    gromacs command when an instance is created.
    """
    # TODO: setup the environment from GMXRC (can use env=DICT in Popen/call)

    command_name = None
    doc_pattern = r'.*?(?P<DOCS>DESCRIPTION.*)'

    def __init__(self,*args,**kwargs):
        """Set up the command with gromacs flags as keyword arguments.::

          cmd = GromacsCommand('v', f=['md1.xtc','md2.xtc'], o='processed.xtc', t=200, ...)

        Gromacs boolean switches (such as ``-v``) are given as python positional
        arguments (``'v'``) or as keyword argument (``v=True``); note the quotes
        in the first case. Negating at boolean switch can be done with
        ``'-nov'``, ``nov=True`` or ``v=False``. Any Gromacs options that take
        parameters are handled as keyword arguments. If an option takes multiple
        arguments (such as the multi-file input ``-f file1 file2 ...``) then the
        list of files must be supplied as a python list.

        The command is executed with the run() method or by
        calling it as a function. The two next lines are equivalent::

          cmd(...)
          cmd.run(...)

        When the command is run one can override options that were given at
        initialization or add additional ones.
        """
        self.gmxargs = self._combineargs(*args, **kwargs)
        self.__doc__ = self.gmxdoc

    def run(self,*args,**kwargs):
        """Run the command; kwargs are added or replace the ones given to the constructor."""
        gmxargs = self.gmxargs.copy()
        gmxargs.update(self._combineargs(*args,**kwargs))
        return self._run_command(**gmxargs)

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
                arglist.append('-no'+flag[1:])  # gromacs switches booleans by prefixing 'no'
            else:
                try:
                    arglist.extend([flag] + value) # option with value list
                except TypeError:
                    arglist.extend([flag, value])  # option with single value
        return map(str, arglist)  # all arguments MUST be strings 

    def _run_command(self,*args,**kwargs):
        """Run command with the given arguments.

           rc,stdout,stderr = c._run_command(*args, input=None, **kwargs)
           
        All positional parameters *args and all gromacs **kwargs are passed on
        to the Gromacs command. input and output keywords allow communication
        with the process.
        
        :Arguments:
        input            string or sequence to be fed to the process' standard input;
                         elements of a sequence are concatenated with
                         newlines, including a railing one    [None]
        stdin            None or automatically set to PIPE if input given [None]
        stdout           filehandle to write to, eg None to see output on screen;
                         PIPE returns the output as a string in the stdout parameter [PIPE]
        stderr           filehandle to write to; STDOUT merges standard error with
                         the standard out stream. PIPE returns the output
                         as a string in the stderr parameter [STDOUT]

        All other kwargs are passed on to the Gromacs tool.
     
        :Returns:
        The shell return code rc of the command is always returned. Depending
        on the value of output, various strings are filled with output from the
        command.

        :Notes:

        At the moment, the processes stdout and stderr are merged        

        STDOUT and PIPE are objects provided by the subprocess module. Any
        python stream can be provided and manipulated. This allows for chaining
        of commands. Use ::
          from subprocess import PIPE, STDOUT
        when requiring the special streams.

        In order to build pipes one must use the special Popen object returned
        by the Popen() method (qv).

        """
        p = self.Popen(*args, **kwargs)
        out, err = p.communicate()       # special Popen knos input!
        rc = p.returncode
        if rc != 0:
            warnings.warn("Command %(command)r failed with return code %(returncode)d" % vars(p),
                          category=RuntimeWarning)
        return rc, out, err

    def Popen(self, *args,**kwargs):
        """Returns a special Popen instance with pre-set input for communicate()."""
        stdin = kwargs.pop('stdin', None)
        stderr = kwargs.pop('stderr', STDOUT)
        stdout = kwargs.pop('stdout', PIPE)     # set to None for screen output
        input = kwargs.pop('input', None)
        if input:
            stdin = PIPE
            if not type(input) is str:
                try:
                    input = '\n'.join(map(str, input)) + '\n'
                except TypeError:
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
        rc,docs,nothing = self._run_command('h', stdout=PIPE)
        m = re.match(self.doc_pattern, docs, re.DOTALL)    # keep from DESCRIPTION onwards
        if m is None:
            return "(No Gromacs documentation available)"
        return m.group('DOCS')
        

    def gmxdoc():
        doc = """Usage for the underlying Gromacs tool (cached)."""
        def fget(self):
            if not (hasattr(self, '__doc_cache') and self.__doc_cache):
                self.__doc_cache = self._get_gmx_docs()
            return self.__doc_cache
        return locals()
    gmxdoc = property(**gmxdoc())

    def help(self,long=True):
        """Print help; same as using ``?`` in ``ipython``."""
        print "\ncommand: %s\n\n" % self.command_name
        print self.__doc__
        if long:
            print "\ncall method: command():\n"
            print self.__call__.__doc__
        
    def __call__(self,*args,**kwargs):
        """Run command with the given arguments.

           rc,stdout,stderr = command(*args, input=None, **kwargs)
           
        All positional parameters *args and all gromacs **kwargs are passed on
        to the Gromacs command. input and output keywords allow communication
        with the process via the python subprocess module.
        
        :Arguments:
        input            string or sequence to be fed to the process' standard input;
                         elements of a sequence are concatenated with
                         newlines, including a railing one    [None]
        stdin            None or automatically set to PIPE if input given [None]
        stdout           filehandle to write to, eg None to see output on screen;
                         PIPE returns the output as a string in the stdout parameter [PIPE]
        stderr           filehandle to write to; STDOUT merges standard error with
                         the standard out stream. PIPE returns the output
                         as a string in the stderr parameter [STDOUT]

        All other kwargs are passed on to the Gromacs tool.
     
        :Returns:
        The shell return code rc of the command is always returned. Depending
        on the value of output, various strings are filled with output from the
        command.

        :Notes:

        By default, the process stdout and stderr are merged.

        In order to chain different commands via pipes one must use the special
        Popen object (see Popen() method of the command) instead of the simple
        call described here and first construct the pipeline explicitly and then
        call the communicate() method of the Popen object.

        STDOUT and PIPE are objects provided by the subprocess module. Any
        python stream can be provided and manipulated. This allows for chaining
        of commands. Use ::
          from subprocess import PIPE, STDOUT
        when requiring the special streams.
        """
        return self.run(*args,**kwargs)




class PopenWithInput(subprocess.Popen):
    """Popen class that knows its input; simply call communicate() later."""
    def __init__(self,*args,**kwargs):
        self.input = kwargs.pop('input',None)
        self.command = args[0]
        super(PopenWithInput,self).__init__(*args,**kwargs)
    def communicate(self):
        return super(PopenWithInput,self).communicate(self.input)


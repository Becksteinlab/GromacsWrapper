# $Id$
"""Core functionality for the Gromacs python shell."""

import subprocess
import warnings

class GromacsCommand(object):
    """Base class for wrapping a g_* command.
    
    Limitations: User must have sourced GMXRC so that the python
    script can inherit the environment and find the gromacs programs.
    """
    # TODO: setup the environment from GMXRC (can use env=DICT in Popen/call)

    command_name = None
    doc_pattern = r'.*?(?P<DOCS>DESCRIPTION.*)'

    def __init__(self,*args,**kwargs):
        """Run the command with gromacs flags as keyword arguments.
        
        result = GromacsCommand('v', f='md.xtc', o='processed.xtc', t=200, ...)
        """
        self.gmxargs = self._args2kwargs(*args)  # switches are treated as kwargs, too
        self.gmxargs.update(kwargs)              # proper options with parameters
        self.__doc__ = self.gmxdoc

    def run(self,*args,**kwargs):
        """Run the command; kwargs are added or replace the ones given to the constructor."""
        return self._run_command(*args,**kwargs)

    def _args2kwargs(self,*args):
        """Turns all args into pseudo-kwargs with value True."""
        return dict([(arg, True) for arg in args])
    
    def _build_arg_list(self,**kwargs):
        """Build list of arguments from the dict; keys must be valid  gromacs flags."""
        arglist = []
        for flag,value in kwargs.items():
            # XXX: check flag against allowed values
            flag = str(flag)
            if not flag.startswith('-'):
                flag = '-' + flag
            # XXX: could treat value == False, too --- useful?
            if value is True:
                arglist.append(flag)          # simple command line flag
            else:
                arglist.extend([flag, value]) # option with value
        return arglist

    def _run_command(self,*args,**kwargs):
        # TODO: flexible in/ouput: 
        # - It should be possible to easily feed input (eg for make_ndx).
        # - Output should be made avialable, maybe even in a
        #   structured manner, depending on the class of the tool.        
        newargs = self._args2kwargs(*args))
        newargs.update(kwargs)        
        gmxkwargs = self.gmxargs.copy()
        gmxargs.update(newargs)
        cmd = [self.command_name] + self._build_arg_list(**gmxargs)
        retcode = subprocess.call(cmd)
        if retcode != 0:
            warnings.warn("Command %(cmd)r failed with return code %(retcode)d" % vars())
        return retcode != 0

    def _get_gmx_docs(self):
        """Extract standard gromacs doc by running the program and chopping the header."""
        cmd = [self.command_name, '-h']
        p = subprocess.Popen(cmd,stderr=subprocess.STDOUT,stdout=subprocess.PIPE)
        docs,nothing = p.communicate()
        m = re.match(self.doc_pattern, docs, re.DOTALL)    # keep from DESCRIPTION onwards
        return m.group('DOCS')

    def gmxdoc():
        doc = """Usage for the underlying Gromacs tool (cached)."""
        def fget(self):
            if not (hasattr(self, '__doc_cache') and self.__doc_cache):
                self.__doc_cache = self._get_gmx_docs()
            return self.__doc_cache
        return locals()
    gmxdoc = property(**gmxdoc())
        
    def __call__(self,*args,**kwargs):
        return self.run(*args,**kwargs)





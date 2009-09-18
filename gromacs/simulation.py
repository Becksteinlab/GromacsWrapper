# $Id$
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)

"""
:mod:`gromacs.simulation` -- Running simulations
================================================

Helper functions and classes around :class:`gromacs.tools.mdrun`.

.. autoclass:: MDrunner

"""

import core
import subprocess
import logging

# logging
logger = logging.getLogger('gromacs.simulation')


class MDrunner(object):
    """A class to run ``mdrun`` in various ways.
    
    In order to do complicated multiprocessor runs with mpiexec or
    similar you need to derive from this class and override

    - :attr:`MDrunner.mdrun` with the path to the ``mdrun`` executable
    - :attr:`MDrunner.mpiexec` with the path to the MPI launcher
    - :meth:`MDrunner.mpicommand` with a function that returns the mpi command as a list
    """

    #: path to the mdrun executable (or the name when on :envvar:`PATH`)
    mdrun = "mdrun"
    #: path to the MPI launcher (e.g. ``mpiexec``)
    mpiexec = None

    def __init__(self, **kwargs):
        """Set up a simple run with ``mdrun``.

        Use the *kwargs* in order to set :class:`gromacs.tools.mdrun`
        options such as *deffnm* = "md" etc.
        """
        # use a GromacsCommand class for handling arguments
        cls = type('MDRUN', (core.GromacsCommand,), 
                   {'command_name': self.mdrun,
                    '__doc__': "MDRUN command %r" % self.mdrun})

        kwargs['failure'] = 'raise'    # failure mode of class
        self.MDRUN = cls(**kwargs)  # might fail for mpi binaries? .. -h?
    
    def commandline(self, **mpiargs):
        """Returns simple command line to invoke mdrun.
        
        Only allows primitive mpi at the moment:
           *mpiexec* -n *ncores* *mdrun* *mdrun-args*
        """
        cmd = self.MDRUN.commandline()
        if self.mpiexec:
            cmd = self.mpicommand(**mpiargs) + cmd
        return cmd

    def mpicommand(self, *args, **kwargs):
        """Return a list of the mpi command portion of the commandline.

        (This is a primitive example for OpenMP. Override it for more
        complicated cases.)
        """
        # example implementation
        ncores = kwargs.pop('ncores', 8)
        return [self.mpiexec, '-n', str(ncores)]

    def run(self, **mpiargs):
        cmd = self.commandline(**mpiargs)
        
        logger.info(" ".join(cmd))
        try:
            rc = subprocess.call(cmd)
        except:
            logger.fatal("Failed MD run for unknown reasons.")
            raise
        if rc == 0:
            logger.info("MDrun completed ok, returncode = %d" % rc)
        else:
            logger.critical("Failure in MDrun, returncode = %d" % rc)

        return rc
        
class MDrunnerOpenMP(MDrunner):
    mdrun = "mdrun_openmp"
    mpiexec = "mpiexec"

class MDrunnerOpenMP64(MDrunner):
    mdrun = "mdrun_openmp64"
    mpiexec = "mpiexec"



def check_mdrun_success(logfile):
    """Check if ``mdrun`` finished successfully.

    Analyses the output from ``mdrun`` in *logfile*. Right now we are
    simply looking for the line "Finished mdrun on node" in the last 1kb of
    the file. (The file must be seeakable.)

    :Arguments:
      logfile : filename
         Logfile produced by ``mdrun``.
    :Returns: boolean (``True`` if all ok, ``False`` otherwise)
    """
    status = False
    log = open(logfile)
    try:
        log.seek(-1024L, 2)
        for line in log:
            if line.startswith("Finished mdrun on node"):
                status = True
                break
    finally:
        log.close()

    return status
    

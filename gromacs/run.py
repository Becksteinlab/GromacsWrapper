# GromacsWrapper: run.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
:mod:`gromacs.run` -- Running simulations
=========================================

Helper functions and classes around :class:`gromacs.tools.Mdrun`.

.. autoclass:: MDrunner
   :members:
.. autoclass:: MDrunnerOpenMP
.. autoclass:: MDrunnerOpenMP64
.. autoclass:: MDrunnerMpich2Smpd

.. autofunction:: check_mdrun_success

"""
from __future__ import absolute_import, with_statement
__docformat__ = "restructuredtext en"

import subprocess
import os.path

# logging
import logging
logger = logging.getLogger('gromacs.run')


# gromacs modules
from . import core
from . import utilities

class MDrunner(utilities.FileUtils):
    """A class to manage running :program:`mdrun` in various ways.

    In order to do complicated multiprocessor runs with mpiexec or
    similar you need to derive from this class and override

    - :attr:`MDrunner.mdrun` with the path to the ``mdrun`` executable
    - :attr:`MDrunner.mpiexec` with the path to the MPI launcher
    - :meth:`MDrunner.mpicommand` with a function that returns the mpi command as a list

    In addition there are two methods named :meth:`prehook` and
    :meth:`posthook` that are called right before and after the
    process is started. If they are overriden appropriately then they
    can be used to set up a mpi environment.

    The :meth:`run` method can take arguments for the
    :program:`mpiexec` launcher but it can also be used to supersede
    the arguments for :program:`mdrun`.

    .. Note:: Changing :program:`mdrun` arguments permanently changes the default arguments 
              for this instance of :class:`MDrunner`. (This is arguably a bug.)
    """

    #: path to the :program:`mdrun` executable (or the name if it can be found on :envvar:`PATH`)
    mdrun = "mdrun"
    #: path to the MPI launcher (e.g. :program:`mpiexec`)
    mpiexec = None

    def __init__(self, dirname=os.path.curdir, **kwargs):
        """Set up a simple run with ``mdrun``.

        :Keywords:
           *dirname*
               Change to this directory before launching the job. Input
               files must be supplied relative to this directory.
           *keywords*
               All other keword arguments are used to construct the
               :class:`~gromacs.tools.mdrun` commandline. Note that only
               keyword arguments are allowed.

        """
        # run MD in this directory (input files must be relative to this dir!)
        self.dirname = dirname

        # use a GromacsCommand class for handling arguments
        cls = type('MDRUN', (core.GromacsCommand,),
                   {'command_name': self.mdrun,
                    '__doc__': "MDRUN command %r" % self.mdrun})

        kwargs['failure'] = 'raise'    # failure mode of class
        self.MDRUN = cls(**kwargs)  # might fail for mpi binaries? .. -h?

        # analyze command line to deduce logfile name
        logname = kwargs.get('g', None)    # explicit
        if logname in (True, None):        # implicit
            logname = 'md'             # mdrun default
            deffnm = kwargs.get('deffnm', None)
            if not deffnm is None:
                logname = deffnm
        self.logname = os.path.realpath(
            os.path.join(self.dirname, self.filename(logname, ext='log')))

    def commandline(self, **mpiargs):
        """Returns simple command line to invoke mdrun.

        If :attr:`mpiexec` is set then :meth:`mpicommand` provides the mpi
        launcher command that prefixes the actual ``mdrun`` invocation:

           :attr:`mpiexec` [*mpiargs*]  :attr:`mdrun` [*mdrun-args*]

        The *mdrun-args* are set on initializing the class. Override
        :meth:`mpicommand` to fit your system if the simple default
        OpenMP launcher is not appropriate.
        """
        cmd = self.MDRUN.commandline()
        if self.mpiexec:
            cmd = self.mpicommand(**mpiargs) + cmd
        return cmd

    def mpicommand(self, *args, **kwargs):
        """Return a list of the mpi command portion of the commandline.

        Only allows primitive mpi at the moment:
           *mpiexec* -n *ncores* *mdrun* *mdrun-args*

        (This is a primitive example for OpenMP. Override it for more
        complicated cases.)
        """
        if self.mpiexec is None:
            raise NotImplementedError("Override mpiexec to enable the simple OpenMP launcher")
        # example implementation
        ncores = kwargs.pop('ncores', 8)
        return [self.mpiexec, '-n', str(ncores)]

    def prehook(self, **kwargs):
        """Called directly before launching the process."""
        return

    def posthook(self, **kwargs):
        """Called directly after the process terminated (also if it failed)."""
        return

    def run(self, pre=None, post=None, mdrunargs=None, **mpiargs):
        """Execute the mdrun command (possibly as a MPI command) and run the simulation.

        :Keywords:
          *pre*
             a dictionary containing keyword arguments for the :meth:`prehook`
          *post*
             a dictionary containing keyword arguments for the :meth:`posthook`
          *mdrunargs*
             a dictionary with keyword arguments for :program:`mdrun` which supersede
             **and update** the defaults given to the class constructor
          *mpiargs*
             all other keyword arguments that are processed by :meth:`mpicommand`
        """

        if pre is None:
            pre = {}
        if post is None:
            post = {}
        if mdrunargs is not None:
            try:
                self.MDRUN.gmxargs.update(mdrunargs)
            except (ValueError, TypeError):
                msg = "mdrunargs must be a dict of mdrun options, not {0}".format(mdrunargs)
                logger.error(msg)
                raise

        cmd = self.commandline(**mpiargs)

        with utilities.in_dir(self.dirname, create=False):
           try:
               self.prehook(**pre)
               logger.info(" ".join(cmd))
               rc = subprocess.call(cmd)
           except:
               logger.exception("Failed MD run for unknown reasons.")
               raise
           finally:
               self.posthook(**post)
        if rc == 0:
            logger.info("MDrun completed ok, returncode = %d" % rc)
        else:
            logger.critical("Failure in MDrun, returncode = %d" % rc)
        return rc

    def run_check(self, **kwargs):
        """Run :program:`mdrun` and check if run completed when it finishes.

        This works by looking at the mdrun log file for 'Finished
        mdrun on node'. It is useful to implement robust simulation
        techniques.

        :Arguments:
           *kwargs* are keyword arguments that are passed on to
           :meth:`run` (typically used for mpi things)

        :Returns:
           - ``True`` if run completed successfully
           - ``False`` otherwise
        """
        rc = None   # set to something in case we ever want to look at it later (and bomb in the try block)
        try:
            rc = self.run(**kwargs)
        except:
            logger.exception("run_check: caught exception")
        status = self.check_success()
        if status:
            logger.info("run_check: Hooray! mdrun finished successfully")
        else:
            logger.error("run_check: mdrun failed to complete run")
        return status

    def check_success(self):
        """Check if :program:`mdrun` finished successfully.

        (See :func:`check_mdrun_success` for details)
        """
        return check_mdrun_success(self.logname)

class MDrunnerDoublePrecision(MDrunner):
    """Manage running :program:`mdrun_d`.
    """
    mdrun = "mdrun_d"

class MDrunnerOpenMP(MDrunner):
    """Manage running :program:`mdrun` as an OpenMP_ multiprocessor job.

    .. _OpenMP: http://openmp.org/wp/
    """
    mdrun = "mdrun_openmp"
    mpiexec = "mpiexec"

class MDrunnerOpenMP64(MDrunner):
    """Manage running :program:`mdrun` as an OpenMP_ multiprocessor job (64-bit executable).

    .. _OpenMP: http://openmp.org/wp/
    """
    mdrun = "mdrun_openmp64"
    mpiexec = "mpiexec"

class MDrunnerMpich2Smpd(MDrunner):
    """Manage running :program:`mdrun` as mpich2_ multiprocessor job with the SMPD mechanism.

    .. _mpich2: http://www.mcs.anl.gov/research/projects/mpich2/
    """
    mdrun = "mdrun_mpich2"
    mpiexec = "mpiexec"

    def prehook(self, **kwargs):
        """Launch local smpd."""
        cmd = ['smpd', '-s']
        logger.info("Starting smpd: "+" ".join(cmd))
        rc = subprocess.call(cmd)
        return rc

    def posthook(self, **kwargs):
        """Shut down smpd"""
        cmd = ['smpd', '-shutdown']
        logger.info("Shutting down smpd: "+" ".join(cmd))
        rc = subprocess.call(cmd)
        return rc



def check_mdrun_success(logfile):
    """Check if ``mdrun`` finished successfully.

    Analyses the output from ``mdrun`` in *logfile*. Right now we are
    simply looking for the line "Finished mdrun on node" in the last 1kb of
    the file. (The file must be seeakable.)

    :Arguments:
      *logfile* : filename
         Logfile produced by ``mdrun``.

    :Returns: ``True`` if all ok, ``False`` if not finished, and
              ``None`` if the *logfile* cannot be opened
    """
    status = False
    try:
        log = open(logfile)
    except IOError:
        return None
    try:
        log.seek(-1024L, 2)
        for line in log:
            if line.startswith("Finished mdrun on node"):
                status = True
                break
    finally:
        log.close()

    return status


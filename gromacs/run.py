# -*- coding: utf-8 -*-
# GromacsWrapper: run.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

""":mod:`gromacs.run` -- Running simulations
=========================================

The :mod:`gromacs.run` module contains tools for launching a Gromacs MD
simulation with :program:`gmx mdrun`. The basic tool is the :class:`MDrunner`
class that customizes how :class:`mdrun<gromacs.tools.Mdrun>` is actually
called. It enables setting a driver such as :program:`mpiexec` for launching
MPI-enabled runs. The :ref:`example-mdrunner-mpi` should make clearer what one
needs to do.

Additionally, :ref:`helpers` are provided to check and manage MD runs.


.. _example-mdrunner-mpi:

Example: How to create your own MDrunner with ``mpiexec -n``
------------------------------------------------------------

- Question: How do I change the GromacsWrapper configuration file so that
  ``mdrun`` gets called with an ``mpiexec -n`` prefix?

- Answer: That's not directly supported but if you just want to change how
  ``mdrun`` is launched then you can create a custom :class:`MDrunner` for this
  purpose.

In many cases, you really only need the path to :program:`mpiexec` and then you
can just derive your own class :class:`MDrunnerMPI`::

  import gromacs.run
  class MDrunnerMPI(gromacs.run.MDrunner):
      \"\"\"Manage running :program:`mdrun` as an MPI multiprocessor job.\"\"\"

      mdrun = "gmx_mpi mdrun"
      mpiexec = "/opt/local/bin/mpiexec"


The full path to the MPI runner :program:`mpiexec` (or :program:`mpirun`) is
stored in the class attribute :attr:`MDrunnerMPI.mpiexec`.

This class can then be used as ::

  mdrun_mpi = MDrunnerMPI(s="md.tpr", deffnm="md")
  rc = mdrun_mpi.run(ncores=16)

Our :class:`MDrunnerMPI` only supports running ``mpiexec -n ncores
gmx mdrun ...``, i.e., only the ``-n ncores`` arguments for :program:`mpiexec`
is supported. If you need more functionality then you need write your own
:meth:`MDrunner.mpicommand` method, which you would add to your own
:class:`MDrunnerMPI` class.

The included :class:`MDrunnerOpenMP` could be used instead of our own
:class:`MDrunnerMPI`; the only difference is that multiple names of MPI-enabled
``mdrun`` binaries are stored as a tuple in the attribute
:attr:`MDrunnerOpenMP.mdrun` so that the class works for old Gromacs 4.x and
modern Gromacs ≥ 2016.

If you need to run some code before or after launching you can add it as the
:meth:`MDrunnerMPI.prehook` and :meth:`MDrunnerMPI.posthook` methods as shown
in :class:`MDrunnerMpich2Smpd`.



MDrunner
--------

The :class:`MDrunner` wraps :class:`gromacs.tools.Mdrun` to customize launching
a Gromacs MD simulation from inside the Python interpreter.

.. autoclass:: MDrunner
   :members:

.. autoclass:: MDrunnerDoublePrecision


Example implementations
-----------------------

.. autoclass:: MDrunnerOpenMP
   :members:
.. autoclass:: MDrunnerMpich2Smpd
   :members:

.. _helpers:

Helper functions
----------------

.. autofunction:: check_mdrun_success
.. autofunction:: get_double_or_single_prec_mdrun
.. autofunction:: find_gromacs_command

"""
__docformat__ = "restructuredtext en"

import warnings
import subprocess
import os.path
import errno
import signal

# logging
import logging

logger = logging.getLogger("gromacs.run")


# gromacs modules
import gromacs
from .exceptions import GromacsError, AutoCorrectionWarning
from . import core
from . import utilities


def find_gromacs_command(commands):
    """Return *driver* and *name* of the first command that can be found on :envvar:`PATH`"""

    # We could try executing 'name' or 'driver name' but to keep things lean we
    # just check if the executables can be found and then hope for the best.

    commands = utilities.asiterable(commands)
    for command in commands:
        try:
            driver, name = command.split()
        except ValueError:
            driver, name = None, command

        executable = driver if driver else name
        if utilities.which(executable):
            break
    else:
        raise OSError(
            errno.ENOENT, "No Gromacs executable found in", ", ".join(commands)
        )

    return driver, name


class MDrunner(utilities.FileUtils):
    """A class to manage running :program:`mdrun` in various ways.

    In order to do complicated multiprocessor runs with mpiexec or similar you
    need to derive from this class and override

    - :attr:`MDrunner.mdrun` with the path to the ``mdrun`` executable
    - :attr:`MDrunner.mpiexec` with the path to the MPI launcher
    - :meth:`MDrunner.mpicommand` with a function that returns the mpi command as a list

    In addition there are two methods named :meth:`prehook` and
    :meth:`posthook` that are called right before and after the process is
    started. If they are overriden appropriately then they can be used to set
    up a mpi environment.

    The :meth:`run` method can take arguments for the :program:`mpiexec`
    launcher but it can also be used to supersede the arguments for
    :program:`mdrun`.

    The actual **mdrun** command is set in the class-level attribute
    :attr:`mdrun`. This can be a single string or a sequence (tuple) of
    strings. On instantiation, the first entry in :attr:`mdrun` that can be
    found on the :envvar:`PATH` is chosen (with
    :func:`find_gromacs_command`). For example, ``gmx mdrun`` from
    Gromacs 5.x but just ``mdrun`` for Gromacs 4.6.x. Similarly, alternative
    executables (such as double precision) need to be specified here
    (e.g. ``("mdrun_d", "gmx_d mdrun")``).

    .. Note:: Changing :program:`mdrun` arguments permanently changes the
              default arguments for this instance of :class:`MDrunner`. (This
              is arguably a bug.)

    .. versionchanged:: 0.5.1
       Added detection of bare Gromacs commands (Gromacs 4.6.x) or commands run through
       :program:`gmx` (Gromacs 5.x).

    .. versionchanged:: 0.6.0
       Changed syntax for Gromacs 5.x commands.
    """

    #: Path to the :program:`mdrun` executable (or the name if it can be found on :envvar:`PATH`);
    #: this can be a tuple and then the program names are tried in sequence. For Gromacs 5
    #: prefix with the driver command, e.g., ``gmx mdrun``.
    #:
    #: .. versionadded:: 0.5.1
    mdrun = ("mdrun", "gmx mdrun")
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
        self.driver, self.name = find_gromacs_command(self.mdrun)

        # use a GromacsCommand class for handling arguments
        cls = type(
            "MDRUN",
            (core.GromacsCommand,),
            {
                "command_name": self.name,
                "driver": self.driver,
                "__doc__": "MDRUN command {0} {1}".format(self.driver, self.name),
            },
        )

        kwargs["failure"] = "raise"  # failure mode of class
        self.MDRUN = cls(**kwargs)  # might fail for mpi binaries? .. -h?

        # analyze command line to deduce logfile name
        logname = kwargs.get("g")  # explicit
        if logname in (True, None):  # implicit
            logname = "md"  # mdrun default
            deffnm = kwargs.get("deffnm")
            if deffnm is not None:
                logname = deffnm
        self.logname = os.path.realpath(
            os.path.join(self.dirname, self.filename(logname, ext="log"))
        )
        self.process = None
        self.signal_handled = False
        signal.signal(signal.SIGINT, self.signal_handler)

    def signal_handler(self, signum, frame):
        """Custom signal handler for SIGINT."""
        if self.process is not None:
            try:
                self.process.terminate()  # Attempt to terminate the subprocess
                self.process.wait()  # Wait for the subprocess to terminate
                self.signal_handled = True
            except Exception as e:
                logger.error(f"Error terminating subprocess: {e}")
        raise KeyboardInterrupt  # Re-raise the KeyboardInterrupt to exit the main script

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
            raise NotImplementedError(
                "Override mpiexec to enable the simple OpenMP launcher"
            )
        # example implementation
        ncores = kwargs.pop("ncores", 8)
        return [self.mpiexec, "-n", str(ncores)]

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
                msg = "mdrunargs must be a dict of mdrun options, not {0}".format(
                    mdrunargs
                )
                logger.error(msg)
                raise

        cmd = self.commandline(**mpiargs)

        with utilities.in_dir(self.dirname, create=False):
            try:
                self.prehook(**pre)
                logger.info(" ".join(cmd))
                self.process = subprocess.Popen(cmd)  # Use Popen instead of call
                returncode = self.process.wait()  # Wait for the process to complete
            except KeyboardInterrupt:
                # Handle the keyboard interrupt gracefully
                logger.info("Keyboard Interrupt received, terminating the subprocess.")
                raise
            except:
                logger.exception("Failed MD run for unknown reasons.")
                raise
            finally:
                self.posthook(**post)
        if returncode == 0:
            logger.info("MDrun completed ok, returncode = {0:d}".format(returncode))
        else:
            logger.critical("Failure in MDrun, returncode = {0:d}".format(returncode))
        return returncode

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
        rc = None  # set to something in case we ever want to look at it later (and bomb in the try block)
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
    """Manage running :program:`mdrun_d`."""

    mdrun = ("mdrun_d", "gmx_d mdrun")


class MDrunnerOpenMP(MDrunner):
    """Manage running :program:`mdrun` as an OpenMP_ multiprocessor job.

    .. _OpenMP: http://openmp.org/wp/
    """

    mdrun = ("mdrun_openmp", "gmx_openmp mdrun")
    mpiexec = "mpiexec"


class MDrunnerMpich2Smpd(MDrunner):
    """Manage running :program:`mdrun` as mpich2_ multiprocessor job with the SMPD mechanism.

    .. _mpich2: http://www.mcs.anl.gov/research/projects/mpich2/
    """

    mdrun = "mdrun_mpich2"
    mpiexec = "mpiexec"

    def prehook(self, **kwargs):
        """Launch local smpd."""
        cmd = ["smpd", "-s"]
        logger.info("Starting smpd: " + " ".join(cmd))
        rc = subprocess.call(cmd)
        return rc

    def posthook(self, **kwargs):
        """Shut down smpd"""
        cmd = ["smpd", "-shutdown"]
        logger.info("Shutting down smpd: " + " ".join(cmd))
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
    if not os.path.exists(logfile):
        return None
    with open(logfile, "rb") as log:
        log.seek(-1024, 2)
        for line in log:
            line = line.decode("ASCII")
            if line.startswith("Finished mdrun on"):
                return True
    return False


def get_double_or_single_prec_mdrun():
    """Return double precision ``mdrun`` or fall back to single precision.

    This convenience function tries :func:`gromacs.mdrun_d` first and
    if it cannot run it, falls back to :func:`gromacs.mdrun` (without
    further checking).

    .. versionadded:: 0.5.1
    """
    try:
        gromacs.mdrun_d(h=True, stdout=False, stderr=False)
        logger.debug("using double precision gromacs.mdrun_d")
        return gromacs.mdrun_d
    except (AttributeError, GromacsError, OSError):
        # fall back to mdrun if no double precision binary
        wmsg = (
            "No 'mdrun_d' binary found so trying 'mdrun' instead.\n"
            "(Note that energy minimization runs better with mdrun_d.)"
        )
        logger.warning(wmsg)
        warnings.warn(wmsg, category=AutoCorrectionWarning)
        return gromacs.mdrun

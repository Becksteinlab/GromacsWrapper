# qsub -- utilities for batch submission systems
# Copyright (c) 2010 Oliver Beckstein <orbeckst@gmail.com>
# Made available under GNU Pulic License v3.

"""
:mod:`gromacs.qsub` -- utilities for batch submission systems
=============================================================

The module helps writing submission scripts for various batch submission
queuing systems. The known ones are listed stored as
:class:`~gromacs.qsub.QueuingSystem` instances in
:data:`~gromacs.qsub.queuing_systems`; append new ones to this list.

The working paradigm is that template scripts are provided (see
:data:`gromacs.config.templates`) and only a few place holders are substituted
(using :func:`gromacs.cbook.edit_txt`).

*User-supplied template scripts* can be stored in
:data:`gromacs.config.qscriptdir` (by default ``~/.gromacswrapper/qscripts``)
and they will be picked up before the package-supplied ones.

The :class:`~gromacs.qsub.Manager` handles setup and control of jobs
in a queuing system on a remote system via :program:`ssh`.

At the moment, some of the functions in :mod:`gromacs.setup` use this module
but it is fairly independent and could conceivably be used for a wider range of
projects.


Queuing system templates
------------------------

The queuing system scripts are highly specific and you will need to add
your own.  Templates should be shell scripts. Some parts of the
templates are modified by the
:func:`~gromacs.qsub.generate_submit_scripts` function. The "place
holders" that can be replaced are shown in the table below. Typically,
the place holders are either shell variable assignments or batch
submission system commands. The table shows SGE_ commands but PBS_ and
LoadLeveler_ have similar constructs; e.g. PBS commands start with
``#PBS`` and LoadLeveller uses ``#@`` with its own command keywords).

.. Table:: Substitutions in queuing system templates.

   ===============  ===========  ================ ================= =====================================
   place holder     default      replacement      description       regex
   ===============  ===========  ================ ================= =====================================
   #$ -N            GMX_MD       *sgename*        job name          `/^#.*(-N|job_name)/`
   #$ -l walltime=  00:20:00     *walltime*       max run time      `/^#.*(-l walltime|wall_clock_limit)/`
   #$ -A            BUDGET       *budget*         account           `/^#.*(-A|account_no)/`
   DEFFNM=          md           *deffnm*         default gmx name  `/^ *DEFFNM=/`
   STARTDIR=        .            *startdir*       remote jobdir     `/^ *STARTDIR=/`
   WALL_HOURS=      0.33         *walltime* h     mdrun's -maxh     `/^ *WALL_HOURS=/`
   NPME=                         *npme*           PME nodes         `/^ *NPME=/`
   MDRUN_OPTS=      ""           *mdrun_opts*     more options      `/^ *MDRUN_OPTS=/`
   ===============  ===========  ================ ================= =====================================

Lines with place holders should not have any white space at the
beginning. The regular expression pattern ("regex") is used to find
the lines for the replacement and the literal default values
("default") are replaced. (Exception: any value that follows an equals
sign "=" is replaced, regardless of the default value in the table
*except* for ``MDRUN_OPTS`` where *only "" will be replace*.) Not all
place holders have to occur in a template; for instance, if a queue
has no run time limitation then one would probably not include
*walltime* and *WALL_HOURS* place holders.

The line ``# JOB_ARRAY_PLACEHOLDER`` can be replaced by
:func:`~gromacs.qsub.generate_submit_array` to produce a "job array"
(also known as a "task array") script that runs a large number of
related simulations under the control of a single queuing system
job. The individual array tasks are run from different sub
directories. Only queuing system scripts that are using the
:program:`bash` shell are supported for job arrays at the moment.

A queuing system script *must* have the appropriate suffix to be properly
recognized, as shown in the table below.

.. Table:: Suffices for queuing system templates. Pure shell-scripts are only used to run locally.

   ==============================  ===========  ===========================
   Queuing system                  suffix       notes
   ==============================  ===========  ===========================
   Sun Gridengine                  .sge         Sun's `Sun Gridengine`_
   Portable Batch queuing system   .pbs         OpenPBS_ and `PBS Pro`_
   LoadLeveler                     .ll          IBM's `LoadLeveler`_
   bash script                     .bash, .sh   `Advanced bash scripting`_
   csh script                      .csh         avoid_ csh_
   ==============================  ===========  ===========================

.. _OpenPBS: http://www.mcs.anl.gov/research/projects/openpbs/
.. _PBS: OpenPBS_
.. _PBS Pro: http://www.pbsworks.com/Product.aspx?id=1
.. _Sun Gridengine: http://gridengine.sunsource.net/
.. _SGE: Sun Gridengine_
.. _LoadLeveler: http://www-03.ibm.com/systems/software/loadleveler/index.html
.. _Advanced bash scripting: http://tldp.org/LDP/abs/html/
.. _avoid: http://www.grymoire.com/Unix/CshTop10.txt
.. _csh: http://www.faqs.org/faqs/unix-faq/shell/csh-whynot/


Example queuing system script template for PBS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following script is a usable PBS_ script for a super computer. It
contains almost all of the replacement tokens listed in the table
(indicated by ++++++). ::

   #!/bin/bash
   # File name: ~/.gromacswrapper/qscripts/supercomputer.somewhere.fr_64core.pbs
   #PBS -N GMX_MD
   #       ++++++
   #PBS -j oe
   #PBS -l select=8:ncpus=8:mpiprocs=8
   #PBS -l walltime=00:20:00
   #                ++++++++

   # host: supercomputer.somewhere.fr
   # queuing system: PBS

   # set this to the same value as walltime; mdrun will stop cleanly
   # at 0.99 * WALL_HOURS
   WALL_HOURS=0.33
   #          ++++

   # deffnm line is possibly modified by gromacs.setup
   # (leave it as it is in the template)
   DEFFNM=md
   #      ++

   TPR=${DEFFNM}.tpr
   OUTPUT=${DEFFNM}.out
   PDB=${DEFFNM}.pdb

   MDRUN_OPTS=""
   #          ++

   # If you always want to add additional MDRUN options in this script then
   # you can either do this directly in the mdrun commandline below or by
   # constructs such as the following:
   ## MDRUN_OPTS="-npme 24 $MDRUN_OPTS"

   # JOB_ARRAY_PLACEHOLDER
   #++++++++++++++++++++++   leave the full commented line intact!

   # avoids some failures
   export MPI_GROUP_MAX=1024
   # use hard coded path for time being
   GMXBIN="/opt/software/SGI/gromacs/4.0.3/bin"
   MPIRUN=/usr/pbs/bin/mpiexec
   APPLICATION=$GMXBIN/mdrun_mpi

   $MPIRUN $APPLICATION -stepout 1000 -deffnm ${DEFFNM} -s ${TPR} -c ${PDB} -cpi \
                        $MDRUN_OPTS \
                        -maxh ${WALL_HOURS} > $OUTPUT
   rc=$?

   # dependent jobs will only start if rc == 0
   exit $rc

Save the above script in ``~/.gromacswrapper/qscripts`` under the name
``supercomputer.somewhere.fr_64core.pbs``. This will make the script
immediately usable. For example, in order to set up a production MD run with
:func:`gromacs.setup.MD` for this super computer one would use ::

   gromacs.setup.MD(..., qscripts=['supercomputer.somewhere.fr_64core.pbs', 'local.sh'])

This will generate submission scripts based on
``supercomputer.somewhere.fr_64core.pbs`` and also the default ``local.sh``
that is provided with *GromacsWrapper*.

In order to modify ``MDRUN_OPTS`` one would use the additonal *mdrun_opts*
argument, for instance::

   gromacs.setup.MD(..., qscripts=['supercomputer.somewhere.fr_64core.pbs', 'local.sh'],
                    mdrun_opts="-v -npme 20 -dlb yes -nosum")


Currently there is no good way to specify the number of processors when
creating run scripts. You will need to provide scripts with different numbers
of cores hard coded or set them when submitting the scripts with command line
options to :program:`qsub`.



Classes and functions
---------------------

.. autoclass:: QueuingSystem
   :members:
.. autofunction:: generate_submit_scripts
.. autofunction:: generate_submit_array
.. autofunction:: detect_queuing_system

.. autodata:: queuing_systems

.. SeeAlso:: :mod:`gromacs.manager` for classes to manage jobs remotely.

"""
from __future__ import absolute_import, with_statement

import os, errno
import warnings

from . import config
from . import cbook
from .utilities import asiterable, Timedelta
from .exceptions import AutoCorrectionWarning

import logging
logger = logging.getLogger('gromacs.qsub')

try:
    from os.path import relpath
except ImportError:
    # appeared in python 2.6
    def relpath(path, start=os.path.curdir):
        """Return a relative version of a path (from posixpath 2.6)"""

        if not path:
            raise ValueError("no path specified")

        start_list = os.path.abspath(start).split(os.path.sep)
        path_list = os.path.abspath(path).split(os.path.sep)

        # Work out how much of the filepath is shared by start and path.
        i = len(os.path.commonprefix([start_list, path_list]))

        rel_list = [os.path.pardir] * (len(start_list)-i) + path_list[i:]
        if not rel_list:
            return os.path.curdir
        return os.path.join(*rel_list)


class QueuingSystem(object):
    """Class that represents minimum information about a batch submission system."""

    def __init__(self, name, suffix, qsub_prefix, array_variable=None, array_option=None):
        """Define a queuing system's functionality

        :Arguments:
          *name*
             name of the queuing system, e.g. 'Sun Gridengine'
          *suffix*
             suffix of input files, e.g. 'sge'
          *qsub_prefix*
             prefix string that starts a qsub flag in a script, e.g. '#$'

        :Keywords:
          *array_variable*
             environment variable exported for array jobs, e.g.
             'SGE_TASK_ID'
          *array_option*
             qsub option format string to launch an array (e.g. '-t %d-%d')

        """
        self.name = name
        self.suffix = suffix
        self.qsub_prefix = qsub_prefix
        self.array_variable = array_variable
        self.array_option = array_option

    def flag(self, *args):
        """Return string for qsub flag *args* prefixed with appropriate inscript prefix."""
        return " ".join((self.qsub_prefix,)+args)

    def has_arrays(self):
        """True if known how to do job arrays."""
        return not self.array_variable is None

    def array_flag(self, directories):
        """Return string to embed the array launching option in the script."""
        return self.flag(self.array_option % (1,len(directories)))

    def array(self, directories):
        """Return multiline string for simple array jobs over *directories*.

        .. Warning:: The string is in ``bash`` and hence the template must also
                     be ``bash`` (and *not* ``csh`` or ``sh``).
        """
        if not self.has_arrays():
            raise NotImplementedError('Not known how make array jobs for '
                                      'queuing system %(name)s' % vars(self))
        hrule = '#'+60*'-'
        lines = [
            '',
            hrule,
            '# job array:',
            self.array_flag(directories),
            hrule,
            '# directories for job tasks',
            'declare -a jobdirs']
        for i,dirname in enumerate(asiterable(directories)):
            idx = i+1   # job array indices are 1-based
            lines.append('jobdirs[%(idx)d]=%(dirname)r' % vars())
        lines.extend([
                '# Switch to the current tasks directory:',
                'wdir="${jobdirs[${%(array_variable)s}]}"' % vars(self),
                'cd "$wdir" || { echo "ERROR: failed to enter $wdir."; exit 1; }',
                hrule,
                ''
                ])
        return "\n".join(lines)

    def isMine(self, scriptname):
        """Primitive queuing system detection; only looks at suffix at the moment."""
        suffix = os.path.splitext(scriptname)[1].lower()
        if suffix.startswith('.'):
            suffix = suffix[1:]
        return self.suffix == suffix

    def __repr__(self):
        return "<"+self.name+" QueuingSystem instance>"

#: Pre-defined queuing systems (SGE, PBS). Add your own here.
queuing_systems = [
    QueuingSystem('Sun Gridengine', 'sge', '#$', array_variable='SGE_TASK_ID', array_option='-t %d-%d'),
    QueuingSystem('PBS', 'pbs', '#PBS', array_variable='PBS_ARRAY_INDEX', array_option='-J %d-%d'),
    QueuingSystem('LoadLeveler', 'll', '#@'),  # no idea how to do arrays in LL
    ]

def detect_queuing_system(scriptfile):
    """Return the queuing system for which *scriptfile* was written."""
    for qs in queuing_systems:
        if qs.isMine(scriptfile):
            return qs
    return None

def generate_submit_scripts(templates, prefix=None, deffnm='md', jobname='MD', budget=None,
                            mdrun_opts=None, walltime=1.0, jobarray_string=None, startdir=None,
                            npme=None, **kwargs):
    """Write scripts for queuing systems.


    This sets up queuing system run scripts with a simple search and replace in
    templates. See :func:`gromacs.cbook.edit_txt` for details. Shell scripts
    are made executable.

    :Arguments:
      *templates*
          Template file or list of template files. The "files" can also be names
          or symbolic names for templates in the templates directory. See
          :mod:`gromacs.config` for details and rules for writing templates.
      *prefix*
          Prefix for the final run script filename; by default the filename will be
          the same as the template. [None]
      *dirname*
          Directory in which to place the submit scripts. [.]
      *deffnm*
          Default filename prefix for :program:`mdrun` ``-deffnm`` [md]
      *jobname*
          Name of the job in the queuing system. [MD]
      *budget*
          Which budget to book the runtime on [None]
      *startdir*
          Explicit path on the remote system (for run scripts that need to `cd`
          into this directory at the beginning of execution) [None]
      *mdrun_opts*
          String of additional options for :program:`mdrun`.
      *walltime*
          Maximum runtime of the job in hours. [1]
      *npme*
          number of PME nodes
      *jobarray_string*
          Multi-line string that is spliced in for job array functionality
          (see :func:`gromacs.qsub.generate_submit_array`; do not use manually)
      *kwargs*
          all other kwargs are ignored

    :Returns: list of generated run scripts
    """
    if not jobname[0].isalpha():
        jobname = 'MD_'+jobname
        wmsg = "To make the jobname legal it must start with a letter: changed to %r" % jobname
        logger.warn(wmsg)
        warnings.warn(wmsg, category=AutoCorrectionWarning)
    if prefix is None:
        prefix = ""
    if not mdrun_opts is None:
        mdrun_opts = '"'+str(mdrun_opts)+'"'  # TODO: could test if quotes already present

    dirname = kwargs.pop('dirname', os.path.curdir)

    wt = Timedelta(hours=walltime)
    walltime = wt.strftime("%h:%M:%S")
    wall_hours = wt.ashours

    def write_script(template):
        submitscript = os.path.join(dirname, prefix + os.path.basename(template))
        logger.info("Setting up queuing system script %(submitscript)r..." % vars())
        # These substitution rules are documented for the user in the module doc string
        cbook.edit_txt(template,
                       [('^ *DEFFNM=','(?<==)(.*)', deffnm),
                        ('^#.*(-N|job_name)', '((?<=-N\s)|(?<=job_name\s))\s*\w+', jobname),
                        ('^#.*(-A|account_no)', '((?<=-A\s)|(?<=account_no\s))\s*\w+', budget),
                        ('^#.*(-l walltime|wall_clock_limit)', '(?<==)(\d+:\d+:\d+)', walltime),
                        ('^ *WALL_HOURS=', '(?<==)(.*)', wall_hours),
                        ('^ *STARTDIR=', '(?<==)(.*)', startdir),
                        ('^ *NPME=', '(?<==)(.*)', npme),
                        ('^ *MDRUN_OPTS=', '(?<==)("")', mdrun_opts),  # only replace literal ""
                        ('^# JOB_ARRAY_PLACEHOLDER', '^.*$', jobarray_string),
                    ],
                       newname=submitscript)
        ext = os.path.splitext(submitscript)[1]
        if ext in ('.sh', '.csh', '.bash'):
            os.chmod(submitscript, 0755)
        return submitscript

    return [write_script(template) for template in config.get_templates(templates)]


def generate_submit_array(templates, directories, **kwargs):
    """Generate a array job.

    For each ``work_dir`` in *directories*, the array job will
     1. cd into ``work_dir``
     2. run the job as detailed in the template
    It will use all the queuing system directives found in the
    template. If more complicated set ups are required, then this
    function cannot be used.

    :Arguments:
       *templates*
          Basic template for a single job; the job array logic is spliced into
          the position of the line ::
              # JOB_ARRAY_PLACEHOLDER
          The appropriate commands for common queuing systems (Sun Gridengine, PBS)
          are hard coded here. The queuing system is detected from the suffix of
          the template.
       *directories*
          List of directories under *dirname*. One task is set up for each
          directory.
       *dirname*
          The array script will be placed in this directory. The *directories*
          **must** be located under *dirname*.
       *kwargs*
          See :func:`gromacs.setup.generate_submit_script` for details.
    """
    dirname = kwargs.setdefault('dirname', os.path.curdir)
    reldirs = [relpath(p, start=dirname) for p in asiterable(directories)]
    missing = [p for p in (os.path.join(dirname, subdir) for subdir in reldirs)
               if not os.path.exists(p)]
    if len(missing) > 0:
        logger.debug("template=%(template)r: dirname=%(dirname)r reldirs=%(reldirs)r", vars())
        logger.error("Some directories are not accessible from the array script: "
                     "%(missing)r", vars())
    def write_script(template):
        qsystem = detect_queuing_system(template)
        if qsystem is None or not qsystem.has_arrays():
            logger.warning("Not known how to make a job array for %(template)r; skipping...", vars())
            return None
        kwargs['jobarray_string'] = qsystem.array(reldirs)
        return generate_submit_scripts(template, **kwargs)[0]   # returns list of length 1

    # must use config.get_templates() because we need to access the file for detecting
    return [write_script(template) for template in config.get_templates(templates)]


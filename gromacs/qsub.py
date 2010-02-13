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
:data:`gromacs.config.templates` for details on the template scripts) and only
a few place holders are substituted (using :func:`gromacs.cbook.edit_txt`).

At the moment, some of the functions in :mod:`gromacs.setup` use this module
but it is fairly independent and could conceivably be used for a wider range of
projects.

.. autoclass:: QueuingSystem
   :members:
.. autofunction:: generate_submit_scripts
.. autofunction:: generate_submit_array
.. autofunction:: detect_queuing_system

.. autodata:: queuing_systems

"""
import os.path
import warnings

import gromacs.config
import gromacs.cbook
from gromacs.utilities import asiterable, Timedelta
from gromacs import AutoCorrectionWarning

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
                            mdrun_opts=None, walltime=1.0, jobarray_string=None,
                            **kwargs):
    """Write scripts for queuing systems.

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
      *mdrun_opts*
          String of additional options for :program:`mdrun`.
      *walltime*
          Maximum runtime of the job in hours. [1]
      *jobarray_string*
          Multi-line string that is spliced in for job array functionality
          (see :func:`gromacs.setup.generate_submit_array`; do not use manually)
      *kwargs*
          all other kwargs are ignored

    :Returns: list of generated run scripts

    This sets up queuing system run scripts with a simple search and replace in
    templates. See :func:`gromacs.cbook.edit_txt` for details.
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
        # These substitution rules are documented for the user in gromacs.config:
        gromacs.cbook.edit_txt(template,
                               [('^ *DEFFNM=','md',deffnm), 
                                ('^#.*(-N|job_name)', 'GMX_MD', jobname),
                                ('^#.*(-A|account_no)', 'BUDGET', budget),
                                ('^#.*(-l walltime|wall_clock_limit)', '00:20:00', walltime),
                                ('^ *WALL_HOURS=', '0\.33', wall_hours),
                                ('^ *MDRUN_OPTS=', '""', mdrun_opts),
                                ('^# JOB_ARRAY_PLACEHOLDER', '^.*$', jobarray_string), 
                                ],
                               newname=submitscript)
        return submitscript

    return [write_script(template) for template in gromacs.config.get_templates(templates)]


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
        logger.error("Some directories are not accessible from the array script: "
                     "%(missing)r" % vars())
    def write_script(template):
        qsystem = detect_queuing_system(template)
        if qsystem is None or not qsystem.has_arrays():
            logger.warning("Not known how to make a job array for %(template)r; skipping..." % vars())
            return None
        kwargs['jobarray_string'] = qsystem.array(reldirs)
        return generate_submit_scripts(template, **kwargs)[0]   # returns list of length 1

    # must use config.get_templates() because we need to access the file for detecting
    return [write_script(template) for template in gromacs.config.get_templates(templates)]

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

The :class:`~gromacs.qsub.Manager` handles setup and control of jobs
in a queuing system on a remote system via :program:`ssh`.

At the moment, some of the functions in :mod:`gromacs.setup` use this module
but it is fairly independent and could conceivably be used for a wider range of
projects.

.. autoclass:: Manager
   :members:
.. autoclass:: QueuingSystem
   :members:
.. autofunction:: generate_submit_scripts
.. autofunction:: generate_submit_array
.. autofunction:: detect_queuing_system

.. autodata:: queuing_systems

Queuing system Manager
----------------------

The :class:`Manager` class must be customized for each system such as
a cluster or a super computer. It then allows submission and control of
jobs remotely (using ssh).


Queuing system templates
------------------------

The queing system scripts are highly specific and you will need to add
your own.  Templates should be sh-scripts and can contain the
following patterns; these are either shell variable assignments or
batch submission system commands. The table shows SGE commands but PBS
and LoadLeveller have similar constructs; e.g. PBS commands start with
``#PBS`` and LoadLeveller uses ``#@`` with its own comman keywords):

===============  ===========  ================ ================= =====================================
command          default      replacement      description       regex
===============  ===========  ================ ================= =====================================
#$ -N            GMX_MD       *sgename*        job name          /^#.*(-N|job_name)/
#$ -l walltime=  00:20:00     *walltime*       max run time      /^#.*(-l walltime|wall_clock_limit)/
#$ -A            BUDGET       *budget*         account           /^#.*(-A|account_no)/
DEFFNM=          md           *deffnm*         default gmx name  /^DEFFNM=/
WALL_HOURS=      0.33         *walltime* h     mdrun's -maxh     /^WALL_HOURS=/
MDRUN_OPTS=      ""           *mdrun_opts*     add.options       /^MDRUN_OPTS=/
===============  ===========  ================ ================= =====================================

These lines should not have any white space at the beginning. The
regular expression pattern is used to find the lines for the
replacement and the default values are replaced.

The line ``# JOB_ARRAY_PLACEHOLDER`` can be replaced by code to run
multiple jobs (a job array) from different sub directories.

"""

import os, errno
import warnings
from subprocess import call, Popen, PIPE
import shutil
import re
import glob

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


class Manager(object):
    """ Base class to launch simulations remotely on computers with queuing systems.

    Basically, ssh into machine and run job.

    Derive a class from :class:`Manager` and override the attributes
     - :attr:`Manager._hostname`
     - :attr:`Manager._scratchdir`
     - :attr:`Manager._qscript`
     - :attr:`Manager._walltime` (if there is a limit to the run time
       of a job)
    and implement a specialized :meth:`Manager.qsub` method if needed.
    """

    # override
    #: hostname of the super computer (**required**)
    _hostname = None
    #: scratch dir on hostname (**required**)    
    _scratchdir = None
    #: name of the template submission script appropriate for the
    #: queuing system on :attr:`Manager._hostname`; can be a path to a
    #: local file or a key for :data:`gromacs.config.templates` (**required**)
    _qscript = None
    #: maximum run time of script in hours; the queuing system script
    #: attr:`Manager._qscript` is supposed to stop :program:`mdrun`
    #: after 99% of this time via the ``-maxh`` option. A value of
    #: ``None`` or ``inf`` indicates to limit.
    _walltime = None

    #: used by :meth:`Manager.get_status`
    log_RE = re.compile(r"""
                Run\stime\sexceeded\s+(?P<exceeded>.*)\s+hours,\swill\sterminate\sthe\srun
                                                                # another part to come
                | Performance:\s*(?P<performance>[\s\d.]+)\n    # performance line  (split later)
                | (?P<completed>Finished\smdrun\son\snode)      # this (part of a) run completed
                """, re.VERBOSE)


    def __init__(self, dirname=os.path.curdir, **kwargs):
        """Set up the manager.

        ssh must be set up (via ~/.ssh/config) to allow access via a
        commandline such as ::

            ssh <hostname> <command> ...

        Typically you want something such as

           host <hostname>
                hostname <hostname>.fqdn.org
                user     <remote_user>

        in ~/.ssh/config and also set up public-key authentication in
        order to avoid typing your password all the time.

        :Arguments:
          *statedir*
              directory component under the remote scratch dir (should
              be different for different jobs)  [basename(CWD)]
          *prefix*
              identifier for job names [MD]
        """
        self.logger = logging.getLogger('gromacs.qsub.Manager')

        # get variables into instance (so that vars(self) works...)
        self.hostname = self._assertnotempty(self._hostname, '_hostname')
        self.scratchdir = self._assertnotempty(self._scratchdir, '_scratchdir')
        self.qscript = self._assertnotempty(self._qscript, '_qscript')
        self.walltime = self._walltime
        self.performance = None   # ns/d, updated with get_status()

        statedir = kwargs.pop('statedir', os.path.basename(os.path.realpath(os.path.curdir)))
        self.wdir = os.path.join(self.scratchdir, statedir)
        self.prefix = kwargs.pop('prefix', 'MD')        # for future use/examples
        self.uri = self.hostname.strip()+":"+self.wdir

        self.logger.info("Setting up a manager from %r.", statedir)

        # test connection and make directory where we run things on the remote host
        rc = call(['ssh', self.hostname, 'mkdir' , '-p', self.wdir])
        if rc == 0:
            self.logger.info("All good: can access %(uri)s" % vars(self))
        else:
            self.logger.error("Problem with ssh and path %(uri)s" % vars(self))

        super(Manager, self).__init__(**kwargs)

    def _assertnotempty(self, value, name):
        """Simple sanity check."""
        if value is None or value == '':
            raise AssertionError("Class %s must have class variable %s set"
                                 % (self.__class__.__name__, name))
        return value

    def remotepath(self, *args):
        """Directory on the remote machine."""
        return os.path.join(self.wdir,*args)

    get_dir = remotepath

    def remoteuri(self, *args):
        """URI of the directory on the remote machine."""
        return os.path.join(self.uri,*args)

    def put(self, dirname):
        """scp dirname to host
        :Arguments: dirname to be transferred
        :Returns: return code from scp
        """
        self.logger.info("Copying %r to %r" % (dirname, self.uri))
	return call(["scp", "-r", dirname, self.uri])

    def putfile(self, filename, dirname):
        """scp *filename* to host in *dirname*
        :Arguments: filename and dirname to be transferred to
        :Returns: return code from scp
        """
        destdir = self.remoteuri(dirname)
        self.logger.info("Copying %(filename)r to %(destdir)r" % vars())
	return call(["scp", filename,  destdir])

    def get(self, dirname, checkfile=None, targetdir=os.path.curdir):
        """``scp -r`` *dirname* from host into *targetdir*

        :Arguments:
        - *dirname*: dir to download
        - *checkfile*: raise OSError/ENOENT if *targetdir/dirname/checkfile* was not found
        - *targetdir*: put *dirname* into this directory
        :Returns: return code from scp
        """
        self.logger.info("Copying %r from %r" % (dirname, self.uri))
        rc = call(["scp", "-r", self.remoteuri(dirname), targetdir])
        #rc = call(["rsync", "-e","ssh","-avP", os.path.join(self.uri,dirname), targetdir])
        if not checkfile is None:
            if not os.path.exists(os.path.join(targetdir, dirname, checkfile)):
                self.logger.error("Failed to get %r from %s", checkfile, self.hostname)
                raise OSError(errno.ENOENT, checkfile, 
                              "Failed to download file from %(hostname)r" % vars(self))
        return rc

    def local_get(self, dirname, checkfile, cattrajectories=True, cleanup=False):
        """Find *checkfile* locally if possible.

        If *checkfile* is not found in *dirname* then it is transferred from the
        remote host.

        If needed, the trajectories are concatenated using :meth:`Manager.cat`.

        :Returns: local path of *checkfile*
        """
        checkpath = os.path.join(dirname, checkfile)
        if not os.path.exists(checkpath):
            self.get(dirname)                # try downloading
            if cattrajectories and not os.path.exists(checkpath):
                # try cating everything first (guess prefix...)
                prefix = os.path.splitext(os.path.basename(checkfile))[0]
                self.cat(dirname, prefix=prefix, cleanup=cleanup)
            if not os.path.exists(checkpath):
                self.logger.error("Failed to get %r from %s", checkfile, self.hostname)
                raise OSError(errno.ENOENT, checkfile, 
                              "Failed to download file from %(hostname)r" % vars(self))
        return checkpath

    def cat(self, dirname, prefix='md', cleanup=True):
        """Concatenate parts of a run in *dirname*.

        Always uses :func:`gromacs.cbook.cat(resolve_multi='guess')`.

        .. Note:: The default is to immediately delete the original files.

        :Keywords:
           *dirname*
              directory to work in
           *prefix*
              prefix (deffnm) of the files [md]
           *cleanup* : boolean
              if ``True``, remove all used files [True]
        """
        gromacs.cbook.cat(prefix, dirname=dirname, resolve_multi='guess')
        # cleanup/get stuff back
        full_dir = os.path.join(dirname, 'full')  # default of cat
        complete_files = os.path.join(full_dir, '*.*')
        self.logger.info("Manager.cat(): recoverning cated files from %r", full_dir)
        for f in glob.glob(complete_files):
            self.logger.debug("Manager.cat(): mv %s %s", f, dirname)
            shutil.move(f, dirname)
        shutil.rmtree(full_dir)
        if cleanup:
            partsdir = os.path.join(dirname, 'parts')  # default of cat
            self.logger.info("Manager.cat(): Removing cat dir %r", partsdir)
            shutil.rmtree(partsdir)
    
    def qsub(self, dirname, **kwargs):
        """Submit remotely on host.

        This is the most simple case: it just runs the command

           cd <remotedir>
           qsub <qscript>
        """

        remotedir = self.remotepath(dirname)
        rc = call(['ssh', self.hostname, 'cd %s && qsub %s' % (remotedir, self.qscript)])
        if rc == 0:
            self.logger.info("Submitted job on %(hostname)s." % vars(self))
        else:
            self.logger.error("Failed running job on on %s in %r.",
                           self.hostname, remotedir)
        return rc == 0

    def get_status(self, dirname, logfilename='md*.log', silent=False):
        """Check status of remote job by looking into the logfile.

        Report on the status of the job and extracts the performance in ns/d if
        available (which is saved in :attr:`Manager.performance`).

        :Arguments:
          - *dirname*
          - *logfilename* can be a shell glob pattern [md*.log]
          - *silent* = True/False; True suppreses log.info messages

        :Returns: ``True`` is job is done, ``False`` if still running
                  ``None`` if no log file found to look at

        .. Note:: Also returns ``False`` if the connection failed.

        .. Warning:; This is an important but **fragile** method. It
                     needs to be improved to be more robust.
        """

        remotefile = os.path.join(self.wdir, dirname, logfilename)

        def loginfo(*args, **kwargs):
            if not silent:
                self.logger.info(*args, **kwargs)
        if not silent:
            self.logger.debug("Checking status of %s:%s", self.hostname, remotefile)

        # need to check if file exists to avoid infinite hangs
        sshcmd = """files=$(ls %(remotefile)s); """ \
            """test -n "$files" && tail -n 500 $(echo $files | tr ' ' '\n' | sort | tail -n 1) """\
            """|| exit 255""" % vars()
        p = Popen(['ssh', self.hostname, sshcmd], 
                  stdout=PIPE, stderr=PIPE, universal_newlines=True)
        out, err = p.communicate()
        rc = p.returncode
        
        status = {'exceeded': False, 'completed': False, 'started': False}
        performance = None
        if rc == 0:
            status['started'] = True
            for m in re.finditer(self.log_RE, out):
                if m.group('completed'):
                    status['completed'] = True
                elif m.group('exceeded'):
                    status['exceeded'] = True
                elif m.group('performance'):
                    performance = dict(zip(['Mnbf/s', 'GFlops', 'ns/day', 'hour/ns'], 
                                           map(float, m.group('performance').split())))
        elif rc == 255:
            loginfo("No output file (yet) for job on %(hostname)s." % vars(self))
            if err:
                self.logger.error("remote: %r", err)
        else:
            self.logger.debug("get_status(): got return code %r, not sure what it means", rc)

        isDone = False
        if status['exceeded']:
            loginfo("Job on %(hostname)s is RUNNING but waiting for next part to run." % vars(self))
        elif status['completed']:  # and not exceeded
            isDone = True
            loginfo("Job on %(hostname)s is DONE." % vars(self))
        elif not status['started']:
            loginfo("Job on %(hostname)s is WAITING in the queue." % vars(self))
        else:
            loginfo("Job on %(hostname)s is still RUNNING." % vars(self))
            if err:
                self.logger.error("remote: %r", err)
            lines = out.split('\n').__iter__()  
            values = ['NAN', 'NAN', 'NAN']   # set a stupid default in case we don't find any time step
            for line in lines:
                if re.match('\s*Step\s+Time\s+Lambda', line):
                    try:
                        values = lines.next().split()  # typically three values
                    except StopIteration:
                        pass                         # if we're unlucky and Step...is last line
            loginfo("Last time step %f ns at step %d.", float(values[1])/1000, float(values[0]))

        if performance:
            self.performance = performance['ns/day']    # used for calculating ndependent()
            loginfo("Performance: %(ns/day)g  ns/day", performance)

        return isDone

    job_done = get_status
    qstat = get_status

    def ndependent(self, runtime,  performance=None, walltime=None):
        """Calculate how many dependent (chained) jobs are required.

        Uses *performance* in ns/d (gathered from :meth:`get_status`) and job max
        *walltime* (in hours) from the class unless provided as keywords.

           n = ceil(runtime/(performance*0.99*walltime)

        :Keywords:
           *runtime*
               length of run in ns
           *performance*
               ns/d with the given setup
           *walltime*
               maximum run length of the script (using 99% of it), in h

        :Returns: *n*  or 1 if walltime is unlimited
        """
        import math
        perf = performance or self.performance    # in ns/d
        wt = walltime or self.walltime            # max runtime of job in h (None = inf)
        
        if wt is None or wt == float('inf'):
            return 1

        if perf is None:
            raise ValueError("No performance data available. Run get_status()?")

        return int(math.ceil(runtime/(perf*0.99*wt/24.)))
        
    def waitfor(self, dirname, **kwargs):
        """Wait until the job associated with *dirname* is done.

        Super-primitive, uses a simple while ... sleep for *seconds* delay

        :Arguments:
          *dirname*
              look for log files under the remote dir corresponding to *dirname*
          *seconds*
              delay in *seconds* during re-polling
        """
        import sys
        import time
        delta_seconds = kwargs.pop('seconds', 120)
        kwargs.setdefault('silent', True)
        totseconds = 0
        while not self.job_done(dirname, **kwargs):
            sys.stderr.write("%4d min   ... %s still running\r" % (totseconds/60, dirname))
            time.sleep(delta_seconds)
            totseconds += delta_seconds
        sys.stderr.write('\n')

    #------------------------------------------------------------
    # example implementations for various stages
    #------------------------------------------------------------
    def setup_posres(self, **kwargs):
        """Set up position restraints run and transfer to host.

        kwargs are passed to :func:`gromacs.setup.MD_restrained`
        
        """
        
        dirname = 'MD_POSRES'
        struct = self.local_get('em','em.pdb')
        gromacs.setup.MD_restrained(dirname=dirname, struct=struct,
                                    qscript=self.qscript, qname=self.prefix+'pr',
                                    **kwargs) 
        self.put(dirname)
        self.logger.info("Run %s on %s in %s/%s" % (dirname, self.hostname, self.uri, dirname))
        self.logger.info(">> qsub('%s')", dirname)
        return dirname

    def setup_MD(self, jobnumber, struct=os.path.join('MD_POSRES', 'md.pdb'), **kwargs):
        """Set up production and transfer to host.
        
        :Arguments:
          - *jobnumber*: 1,2 ...
          - *struct* is the starting structure (default from POSRES
            run but that is just a guess);
          - kwargs are passed to :func:`gromacs.setup.MD`
        """
        kwargs.setdefault('runtime', 1e4)

        jobid_s = '%(jobnumber)03d' % vars()
        dirname = 'MD_'+jobid_s
        structure = self.local_get(os.path.dirname(struct), os.path.basename(struct))

        gromacs.setup.MD(dirname=dirname, struct=structure, qscript=self.qscript, 
                         qname=self.prefix+jobid_s,
                         **kwargs) 
        self.put(dirname)
        self.logger.info("Run %s on %s in %s/%s" % (dirname, self.hostname, self.uri, dirname))
        self.logger.info("Or use %s.qsub(%r)" % (self.__class__.__name__, dirname))

        return dirname

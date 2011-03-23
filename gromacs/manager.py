"""
Managing jobs remotely
======================

.. Warning:: Experimental code, use at your own risk.

Configuration file
------------------

Example::

   [DEFAULT]
   name = leviathan

   [local]   
   topdir = ~

   [remote]
   hostname = leviathan.petagrid.org
   scratchdir = /scratch/username/Projects

   [queuing_system]
   name = PBS
   qscript = leviathan.pbs
   walltime = 24.0

All entries except *walltime* are required; *walltime* can be omitted or set to ``None``.

The remote directory name is constructed in the following way:

1. *topdir* is stripped from the local working directory to give WDIR
2. *scratchdir*/WDIR is the directory on the remote system


See :class:`Manager` for how these values are used.

:func:`get_manager` creates a :class:`Manager` from a configuration file.

Queuing system Manager
----------------------

The :class:`Manager` class must be customized for each system such as
a cluster or a super computer. It then allows submission and control of
jobs remotely (using ssh_).

.. autoclass:: Manager
   :members:
   :exclude-members: job_done, qstat

   .. autoattribute:: _hostname
   .. autoattribute:: _scratchdir
   .. autoattribute:: _qscript
   .. autoattribute:: _walltime
   .. method:: job_done

               alias for :meth:`get_status`

   .. method:: qstat

               alias for :meth:`get_status`


.. _ssh: http://www.openssh.com/
.. _~/.ssh/config: http://linux.die.net/man/5/ssh_config

"""
from __future__ import with_statement

import os
import errno
from subprocess import call, Popen, PIPE
import shutil
import fnmatch
import re
import glob
from ConfigParser import SafeConfigParser

from gromacs import MissingDataError
import gromacs.config
import gromacs.utilities
import gromacs.cbook

import warnings
import logging
logger = logging.getLogger("gromacs.manager")

class ManagerConfigParser(SafeConfigParser):
    def getpath(self, section, option):
        """Return option as an expanded path."""
        return os.path.expanduser(os.path.expandvars(self.get(section, option)))
    def getfloat(self, section, option):
        """Return as :func:`float` or ``None``."""
        try:
            return float(option)
        except ValueError:
            return None

def find_manager_config(name):
    """Find a configuration file for manager *name*."""
    found = list(gromacs.utilities.find_files(gromacs.config.managerdir, name+".cfg"))
    if len(found) == 0:
        errmsg = "No configuration file found for name %r" % name
        logger.error(errmsg)
        raise MissingDataError(errmsg)
    elif len(found) > 1:
        logger.warn("Multiple configuration files found: %r. Using the first one!", found)
    return found[0]

def get_manager_config(filename):
    """Load manager configuration file from *filename*."""
    logger.info("Loading Manager configuration data from %(filename)r", vars())
    cfg = ManagerConfigParser()
    cfg.add_section("local")
    cfg.set("local", "topdir", os.path.expanduser("~"))
    cfg.add_section("remote")
    cfg.add_section("queuing_system")
    cfg.set("queuing_system", "walltime", "None")
    cfg.set("DEFAULT", "filename", filename)
    cfg.readfp(open(filename))
    return cfg

class Job(dict):
    """Properties of a job."""

class Manager(object):
    """ Base class to launch simulations remotely on computers with queuing systems.

    Basically, ssh into machine and run job.

    Derive a class from :class:`Manager` and override the attributes

     - :attr:`Manager._hostname` (hostname of the machine)
     - :attr:`Manager._scratchdir` (all files and directories will be created under 
       this scratch directory; it must be a path on the remote host)
     - :attr:`Manager._qscript` (the default queuing system script template)
     - :attr:`Manager._walltime` (if there is a limit to the run time
       of a job; in hours)

    and implement a specialized :meth:`Manager.qsub` method if needed.

    ssh_ must be set up (via `~/.ssh/config`_) to allow access via a
    commandline such as ::

        ssh <hostname> <command> ...

    Typically you want something such as ::

       host <hostname>
            hostname <hostname>.fqdn.org
            user     <remote_user>

    in ``~/.ssh/config`` and also set up public-key authentication in
    order to avoid typing your password all the time.
    """

    #: Regular expression used by :meth:`Manager.get_status` to parse
    #: the logfile from :program:`mdrun`.
    log_RE = re.compile(r"""
                Run\stime\sexceeded\s+(?P<exceeded>.*)\s+hours,\swill\sterminate\sthe\srun
                                                                # another part to come
                | Performance:\s*(?P<performance>[\s\d.]+)\n    # performance line  (split later)
                | (?P<completed>Finished\smdrun\son\snode)      # this (part of a) run completed
                """, re.VERBOSE)

    def __init__(self, name, dirname=None, **kwargs):
        """Set up the manager.

        :Arguments:
          *name*
              configuration name (corresponds to a store cfg file)
          *dirname*
              directory component under the remote scratch dir (should
              be different for different jobs); the default is to strip
              *topdir* from the config file from the full path of the 
              current directory
          *prefix*
              identifier for job names [MD]
        """
        self.name = name
        self.logger = logging.getLogger('gromacs.manager.%(name)s' % vars())

        try:
            cfg = get_manager_config(find_manager_config(name))
        except:
            logger.error("Failed to read the configuration for Manager %(name)r.", vars())
            raise
        attribs = {
            'name': cfg.get('DEFAULT', 'name'),
            'topdir':  cfg.getpath('local', 'topdir'),
            'hostname': cfg.get('remote', 'hostname'),
            'scratchdir': cfg.get('remote', 'scratchdir'),
            'queuing_system': cfg.get('queuing_system', 'name'),
            'qscript': cfg.get('queuing_system', 'qscript'),
            'walltime': cfg.getfloat('queuing_system', 'walltime'),
            }
        if attribs['name'] != self.name:
            errmsg = "Sanity check failed: requested name %r does not match the Manager name %r "\
                "that was recorded in the config file %r." % \
                (self.name, attribs['name'], cfg.getpath('DEFAULT', 'filename'))
            logger.fatal(errmsg)
            raise ValueError(errmsg)

        self.__dict__.update(attribs)

        self.performance = None   # ns/d, updated with get_status()

        if dirname is None:
            logger.info("Stripping %(topdir)r from current dirname to generate workdir path.", vars(self))
            dirname = os.path.realpath(os.path.curdir).replace(os.path.normpath(self.topdir)+"/", "")
        self.wdir = os.path.normpath(os.path.join(self.scratchdir, dirname))
        self.prefix = kwargs.pop('prefix', 'MD')        # for future use/examples
        self.uri = self.hostname.strip()+":"+self.wdir

        self.logger.info("Setting up a manager from %r.", dirname)

        # test connection and make directory where we run things on the remote host
        rc = call(['ssh', self.hostname, 'mkdir' , '-p', self.wdir])
        if rc == 0:
            self.logger.info("All good: can access %(uri)s" % vars(self))
        else:
            self.logger.error("Problem with ssh and path %(uri)s" % vars(self))

        super(Manager, self).__init__(**kwargs)


    def remotepath(self, *args):
        """Directory on the remote machine."""
        return os.path.join(self.wdir,*args)

    get_dir = remotepath

    def remoteuri(self, *args):
        """URI of the directory on the remote machine."""
        return os.path.join(self.uri,*args)

    def put(self, dirname):
        """scp dirname to host.

        :Arguments: dirname to be transferred
        :Returns: return code from scp
        """
        self.logger.info("Copying %r to %r" % (dirname, self.uri))
	return call(["scp", "-r", dirname, self.uri])

    def putfile(self, filename, dirname):
        """scp *filename* to host in *dirname*.

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

        Always uses :func:`gromacs.cbook.cat` with *resolve_multi* = 'guess'.

        .. Note:: The default is to immediately delete the original files
                  (*cleanup* = ``True``).

        :Keywords:
           *dirname*
              directory to work in
           *prefix*
              prefix (deffnm) of the files [md]
           *cleanup* : boolean
              if ``True``, remove all used files [``True``]
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
        """Submit job remotely on host.

        This is the most primitive implementation: it just runs the commands ::

           cd remotedir && qsub qscript

        on :attr:`Manager._hostname`. *remotedir* is *dirname* under
        :attr:`Manager._scratchdir` and *qscript* defaults to the queuing system
        script hat was produced from the template :attr:`Manager._qscript`.
        """

        remotedir = kwargs.pop('remotedir', self.remotepath(dirname))
        qscript = kwargs.pop('qscript', os.path.basename(self.qscript))
        rc = call(['ssh', self.hostname, 'cd %s && qsub %s' % (remotedir, qscript)])
        if rc == 0:
            self.logger.info("Submitted job %r on %s.", qscript, self.hostname )
        else:
            self.logger.error("Failed running job %s on %s in %r.",
                              qscript, self.hostname, remotedir)
        return rc == 0

    def get_status(self, dirname, logfilename='md*.log', silent=False):
        """Check status of remote job by looking into the logfile.

        Report on the status of the job and extracts the performance in ns/d if
        available (which is saved in :attr:`Manager.performance`).

        :Arguments:
          - *dirname*
          - *logfilename* can be a shell glob pattern [md*.log]
          - *silent* = True/False; True suppresses log.info messages

        :Returns: ``True`` is job is done, ``False`` if still running
                  ``None`` if no log file found to look at

        .. Note:: Also returns ``False`` if the connection failed.

        .. Warning:: This is an important but somewhat  **fragile** method. It
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

        *kwargs* are passed to :func:`gromacs.setup.MD_restrained`
        
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
        
    

# def get_manager(name):
#     """Factory function that creates a new Manager class, based on a config file.

#     :Arguments:
#        *name*
#           name of the config file, `~/.gromacswrapper/managers/name.cfg`
#     """    
#     import qsub
#     warnings.warn("Old-style, derives from qsub.Manager", DeprecationWarning)
#     cfg = get_manager_config(find_manager_config(name))
#     attribs = {
#         'name': cfg.get('manager', 'name'),
#         '_hostname': cfg.get('manager', 'hostname'),
#         '_scratchdir': cfg.get('manager', 'scratchdir'),
#         'queuing_system': cfg.get('queuing_system', 'name'),
#         '_qscript': cfg.get('queuing_system', 'qscript'),
#         '_walltime': cfg.getfloat('queuing_system', 'walltime'),
#         }
#     return type(name, (qsub.Manager,), attribs)


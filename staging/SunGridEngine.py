# $Id: SunGridEngine.py 2765 2009-01-20 13:02:14Z oliver $
"""
:mod:`staging.SunGridEngine` --- staging class for SunGridEngine
================================================================

Primitive framework for staging jobs in `Sun Grid Engine`_ via a
customized :class:`Job` class.

Example python submission script
--------------------------------

Write the SGE script like this::

   #!/usr/bin/env python
   #$ -N bulk
   #$ -S /usr/bin/python
   #$ -v PYTHONPATH=/home/oliver/Library/python-lib
   #$ -v LD_LIBRARY_PATH=/opt/intel/cmkl/8.0/lib/32:/opt/intel/itc60/slib:/opt/intel/ipp41/ia32_itanium/sharedlib:/opt/intel/ipp41/ia32_itanium/sharedlib/linux32:/opt/intel/fc/9.0/lib:/opt/intel/cc/9.0/lib
   #$ -r n
   #$ -j y
   # The next line is IMPORTANT when you are using the default for Job(startdir=None)
   #$ -cwd

   from staging.SunGridEngine import Job

   job = Job(inputfiles=dict(psf = 'inp/crbp_apo.psf',
                             dcd = 'trj/rmsfit_1opa_salt_ewald_shake_10ang_prod.dcd'),
             outputfiles=dict(dx = '*.dx', pickle = '*.pickle'),
             variables=dict(normalize = True, ...))

   job.stage()
   F = job.filenames  # use F[key] to reference filenames from inputfiles or outputfiles
   V = job.variables  # and V[key] for the variables

   # your python script here...
   print "psf: %(psf)s  dcd: %(dcd)" % F
   print "normalize = %(normalize)s" % V


   job.unstage()
   job.cleanup()   # removes stage dir, careful!

.. _`Sun Grid Engine`: http://gridengine.sunsource.net/

Description of the :class:`Job` class
-------------------------------------

.. autoclass:: Job
   :members:

Helper functions for building job arrays
----------------------------------------

.. autofunction:: getline_from_arraylist
.. autofunction:: get_fields_from_arraylist
.. autofunction:: get_value_from_arraylist
"""

import os
import errno
import shutil
from staging.common import joindicts,pathjoin

# TODO: SGE_job should really inherit from Job so that one could
# derive Jobs classes for different queuing systems.

class SGE_job(object):
    """Specifics for a Sun Gridengine job."""

    def __init__(self,*args,**kwargs):
        """Set up a SGE job.

        If the environment contains JOB_NAME, JOB_ID, and SGE_TASK_ID
        then this is a job submitted through Sun Gridengine (we only
        check for JOB_NAME) and staging proceeds as usual.

        If there is no JOB_NAME this is a 'Local' job and we do not do
        any staging, just providing the same framework but on local
        files.

        If the JOB_NAME='my_identifier' *keyword argument* is given
        then staging proceeds as if this was a regular job; in this
        case one can also supply SGE_TASK_ID.

        Arguments:

        JOB_NAME    force use of this JOB_NAME and force staging even
                    when not submitting throgh Gridengine (no effect when SGE controlled)
        SGE_TASK_ID fake a task id in an SGE array job
        """
        super(SGE_job,self).__init__()
        self.__MODE__ = "init"
        if 'JOB_NAME' in os.environ:  # SGE submission
            self.queuingsystem = 'SGE'
            self.JOB_NAME = os.environ['JOB_NAME'].strip()
            self.JOB_ID = os.environ['JOB_ID'].strip()
            self.TASK_ID = os.environ['SGE_TASK_ID'].strip()
        elif 'JOB_NAME' in kwargs:
            self.queuingsystem = 'Copy'
            self.JOB_NAME = kwargs['JOB_NAME'] or "noname"
            self.JOB_ID = str(os.getpid())   # potentially unsafe, use hash or mktemp?
            self.TASK_ID = str(kwargs.setdefault('SGE_TASK_ID','undefined'))
        else:         # running job locally (shouldn't this be in staging.Local?)
            self.queuingsystem = 'Local'
            self.JOB_NAME = 'Local'
            self.JOB_ID = str(os.getpid())
            self.TASK_ID = str(kwargs.setdefault('SGE_TASK_ID','undefined'))
        self.jobdir_name = self.get_jobdir_name()
        self.hostname = self.get_hostname()
        self.__continue_this_line = False     # for msg(), helper for continuing lines

    def get_jobdir_name(self):
        # create canonical name
        if self.TASK_ID == "undefined":
            jobdir_name=self.JOB_NAME+'.'+str(self.JOB_ID)
        else:
            # part of an array job
            jobdir_name=self.JOB_NAME+'.'+str(self.JOB_ID)+'.'+str(self.TASK_ID)
        return jobdir_name

    def get_hostname(self):
        import socket
        return socket.gethostname()

    def msg(self,string,newline=True):
        # suppress newline here
        if not self.__continue_this_line:
            print "%s(): "%self.__MODE__ + str(string),
        else:
            print str(string),
        # add newline if requested
        if newline:
            print
        # tell next invocation what to do
        self.__continue_this_line = not newline

    def statusmessage(self):
        self.msg("hostname:       %s" % self.hostname)
        self.msg("queuing system: %s" % self.queuingsystem)        
        self.msg("JOB_NAME:       %s" % self.JOB_NAME)
        self.msg("JOB_ID:         %s" % self.JOB_ID)
        self.msg("SGE_TASK_ID:    %s" % self.TASK_ID)
        self.msg("jobdir_name:    %s" % self.jobdir_name)


class Job(SGE_job):
    """The Job class encapsulates the SGE job and allows for clean staging and unstaging.

    Set up the Job::

      job = Job(inputfiles=dict(...),outputfiles=dict(...),variables=dict(...),**kwargs)

    *inputfiles* and *outputfiles* are dictionaries with arbitrary
    keys; each item is a path to a file relative to the startdir
    (which by default is the directory from which the SGE job starts
    --- use the ``#$ -cwd`` flag!). If the files are not relative to the
    start dir then new directories are constructed under the stage
    dir; in this instance it uis important that the user script *only*
    uses the filenames in :attr:`Job.filenames`: These have the proper paths
    of the local (staged) files for the script to operate on.

    With ::

      job.stage()

    inputfiles are copied to the stagedir on the node's scratch
    dir and sub directories are created as necessary; directories
    mentioned as part of the outputfiles are created, too. ::

      job.unstage()

    copies back all files mentioned in output files (again, use
    directories as part of the path as necessary) and create the
    directories in the startdir if needed. For the outputfiles one
    can also use shell-style glob patterns, e.g. ``outfiles =
    {'all_dcd': '*.dcd', 'last_data':'*[5-9].dat'}``

    Sensible defaults are automatically selected for startdir
    (cwd) and stagedir (/scratch/USER/JOB_NAME.JOB_ID).

    If the script is not run through SGE (i.e. the environment
    variable :envvar:`JOB_NAME` is not set) then the script is run without
    staging; this is pretty much equivalent to using ::

      from staging.Local import Job

    (i.e. using the :class:`staging.Local.Job` class).

    :Attributes:

       :attr:`input`
                        inputfiles dict  (relative to startdir or absolute)
       :attr:`output`
                        outputfiles dict (relative to startdir or absolute, can contain globs)
       :attr:`filenames`
                        merged dict of input and output, pointing to *staged* files
       :attr:`variables`
                        variables dict

    :Methods:

       :meth:`stage`
                        setup job on the nodes in stagedir
       :meth:`unstage`
                        retrieve results to startdir
       :meth:`cleanup`
                        remove all files on the node (rm -rf stagedir)
    """

    def __init__(self,*args,**kwargs):
        """Set up SGE job.

        :Arguments:

           inputfiles
                            dict of input files (with relative path to startdir);
                            globs are not supported.
           outputfiles
                            dict of result files or glob patterns (relative to
                            stagedir == relative to startdir)
           variables
                            key/value pairs that can be used in the script as 
                            Job.variables[key]
           startdir
                            path to the directory where the input can be found
                            (must be nfs-mounted on node)
           stagedir
                            local scratch directory on node; all input files are copied
                            there. The default should be ok.

           JOB_NAME
                            unique identifier (only set this if this NOT submitted through
                            the Gridengine queuing system AND if the files should be copied
                            to a scratch disk (i.e. staging proceeds as it would for a
                            SGE-submitted job).)
           SGE_TASK_ID
                           fake a task id (use with JOB_NAME)

        """
        self.__MODE__ = "init"   # current state, for self.msg
        super(Job,self).__init__(*args,**kwargs)
        self.input = kwargs.setdefault('inputfiles',{})
        self.output = kwargs.setdefault('outputfiles',{})
        
        self.variables = kwargs.setdefault('variables',{})
        # where we find input files and copy back results
        self.startdir = self.startdir_name(kwargs.setdefault('startdir',None))
        # local directory on node
        self.stagedir = self.stagedir_name(kwargs.setdefault('stagedir',None))
        # normalized filenames (always under stagedir)
        self.filenames = {k: pathjoin(self.stagedir,path,refdir=self.startdir)
             for k,path in joindicts(self.input,self.output).items()}
        
        self.statusmessage()

    def statusmessage(self):
        super(Job,self).statusmessage()        
        self.msg("startdir:       %s" % self.startdir)
        self.msg("stagedir:       %s" % self.stagedir)

    def startdir_name(self,startdir=None):
        if startdir is None:
            # use canonical setup (relies on -cwd SGE flag)
            startdir=os.path.realpath(os.path.curdir)
        return startdir

    def stagedir_name(self,stagedir=None):
        if self.queuingsystem is 'Local':
            return None
        if stagedir is None:
            # use canonical setup
            stagedir = pathjoin('/scratch',os.environ['USER'],self.jobdir_name)
        return stagedir

    def stage(self):
        """Copy all input files to the scratch directory."""
        self.__MODE__ = "stage"
        if self.queuingsystem is 'Local':
            return
        stagedir = self.stagedir
        try:
            os.makedirs(stagedir)
            self.msg("Created stage dir %(stagedir)s." % locals())
        except os.error,e:
            if e.errno == errno.EEXIST:
                self.msg("WARNING %(stagedir)s already exists." % locals())
            else:
                raise                        
        self._make_all_dirs(stagedir,self.input,refdir=self.startdir)  # copy input and preserve directory structure
        self._make_all_dirs(stagedir,self.output,refdir=self.startdir) # also create directories for the output files
        for key,p in self.input.items():                 # copy input files
            srcpath = pathjoin(self.startdir,p, sanitize=False)   # may be absolute (and ignores startdir!)
            destpath = self.filenames[key]               # ALWAYS under stagedir
            self.msg("item=%(key)s: copying %(srcpath)s" % locals(), newline=False)            
            shutil.copyfile(srcpath,destpath)
            self.msg(" --> %(destpath)s" % locals())
        # finally, change current directory to the stage dir: all further
        # commands can assume that staging has been completed            
        os.chdir(stagedir)
        self.msg("chdir to %(stagedir)s successful." % locals())
               
    def unstage(self):
        """Copy results back. Shell-style glob patterns are allowed."""
        self.__MODE__ = "unstage"
        if self.queuingsystem is 'Local':
            return
        import glob
        self._make_all_dirs(self.startdir,self.output,sanitize=False)  # make result directories, may be absolute!
        for key,p in self.output.items():
            src = self.filenames[key]                                 # always relative to stagedir            
            srcdir = os.path.dirname(p)
            destdir = pathjoin(self.startdir,srcdir, sanitize=False)  # may be absolute
            self.msg("item=%(key)s: looking for %(p)s [=%(src)s]..." % locals())
            for srcpath in glob.glob(src):
                srcname = os.path.basename(srcpath)
                destpath = pathjoin(destdir,srcname, sanitize=False)
                self.msg("item=%(key)s: copying %(srcpath)s" % locals(), newline=False)
                shutil.copyfile(srcpath,destpath)   # silently replaces files !
                self.msg(" --> %(destpath)s" % locals())
		
    def cleanup(self):
        """Remove stage dir"""
        self.__MODE__ = "cleanup"
        os.chdir(self.startdir)
        if self.queuingsystem is 'Local':
            return        
        try:
            shutil.rmtree(self.stagedir)
            self.msg("removed stage dir %s" % self.stagedir)
        except os.error,e:
            if e.errno == errno.ENOENT:
                self.msg("%s does not exist any more" % self.stagedir)
            else:
                raise
 
    def _make_all_dirs(self,topdir,filedict,**kwargs):
        """Create directories under topdir, based on paths in filedict."""
        for key,p in filedict.items():
            srcdir = os.path.dirname(p)
            destdir = pathjoin(topdir,srcdir,**kwargs)
            try:
                os.makedirs(destdir)  # recursive
                self.msg("item=%(key)s: created dir %(destdir)s" % locals())
            except os.error,e:
                if e.errno == errno.EEXIST:
                    pass
                else:
                    raise
            
    def save(self,filename):
        """Save the Job() as a pickled file.

        Restore with ::

           import staging.SunGridengine
           import cPickle
           job = cPickle.load(open(<filename>,'r'))
        """
        import cPickle
        cPickle.dump(self,open(filename,'wb'),cPickle.HIGHEST_PROTOCOL)


def getline_from_arraylist(filename=None,ENVNAME='ARRAYLIST',default="arraylist.txt"):
    """Read a list of values from filename and return the line that corresponds to the current SGE_TASK_ID.

      line = get_line_from_arraylist(filename=None,ENVNAME='ARRAYLIST',default="arraylist.txt")

    fields will be different depending on the value of :envvar:`SGE_TASK_ID`
    (set by SunGridengine).  The lines are simply numbered consecutively.

    :Arguments:
       *filename*
                     name of the arraylist file
       *ENVNAME*
                     try to get filename from environment variable if filename is not set
       *default*
                     if all fails, try this as a default filename

    File format::

      # comment lines are ignored as are whitespace lines
      # only the first column is read; the internal numbering starts at 1
      line1 ...   <---- task id 1
      line2 ...   <---- task id 2
      # more comments, they are NOT counted for the task id
      line3 ...   <---- task id 3      
      ...
      
    Ignores white space lines and lines starting with ``#``. Lines are
    stripped of left and right white space.
    """
    if filename is None:
        filename = os.environ.setdefault(ENVNAME, default)
    values = {}
    ival = 0
    # read in file list as dict, indexed by taskid, hence one-based
    for line in open(filename):
        line = line.strip()
        if len(line) == 0 or line[0] == '#':
            continue
        ival +=1
        values[ival] = line

    # print values
    # this varies from task to task (and raises an exception if this is not an array job)
    TASK_ID = os.environ['SGE_TASK_ID']
    if TASK_ID == "undefined":
        raise RuntimeError("This must be run from a SGE task array job.")
    return values[int(TASK_ID)]

def get_fields_from_arraylist(**kwargs):
    """Read a list of values from filename and return the line that corresponds to the current SGE_TASK_ID.

      get_line_from_arraylist(filename=None,ENVNAME='ARRAYLIST',default="arraylist.txt") -> fields

    fields will be different depending on the value of SGE_TASK_ID (set by SunGridengine).
    The lines are simply numbered consecutively.

    See :func:`getline_from_arraylist` for more details.
    """
    return getline_from_arraylist(**kwargs).split()
    
def get_value_from_arraylist(index=0,**kwargs):
    """Get field[index] of the entry in the array list corresponding to SGE_TASK_ID.

    See :func:`get_fields_from_arraylist` for details.
    """
    return get_fields_from_arraylist(**kwargs)[index]
    


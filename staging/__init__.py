# $Id: __init__.py 2677 2008-12-22 20:58:45Z oliver $

"""The 'staging' module provides a framework to run python scripts
easily through a queuing system that requires copying of files to the
scratch directory on compute nodes.

Load the appropiate submodule at the of the script. Currently
available submodules are:

   SunGridEngine       Sun Gridengine (SGE)
                       activates if JOB_ID is found in environment
                       (otherwise behaves similarly to Local)
   Local               Runs the job without any staging.


Example:

  from staging.SunGridEngine import Job

  job = Job(inputfiles=dict(psf = 'inp/apo.psf',
                          dcd = 'trj/prod.dcd'),
            outputfiles=dict(dx = '*.dx', pickle = '*.pickle'))

  job.stage()     # copy files to staging directory, creating dirs

  # your python script here...
  # access all filenames as job.filenames[<KEY>]

  job.unstage()   # copies files backs and creates dirs as needed
  job.cleanup()   # removes stage dir, careful!
"""

__all__ = [ 'SunGridEngine', 'Local' ]

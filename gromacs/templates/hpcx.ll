#  $Id$
#@ shell = /bin/bash
#
#@ job_name = GMX_MD
#
#@ job_type = parallel
#@ cpus = 64
#@ requirements = ( Feature == "SMT" )
#@ tasks_per_node = 32
#@ node_usage = not_shared
#
#@ network.MPI = csss,shared,US
#@ bulkxfer = yes
#
#@ wall_clock_limit = 00:20:00
#@ account_no = BUDGET
#
#@ output = $(job_name).$(schedd_host).$(jobid).out
#@ error  = $(job_name).$(schedd_host).$(jobid).err
#@ notification = never
#
#@ queue

# Basic LoadLeveller script for HPCx from
# http://www.hpcx.ac.uk/support/documentation/UserGuide/HPCxuser/Batch_Processing.html#SECTION00082100000000000000
# - SMT is enabled as it gives a 1/3 speedup for free (tested on a 65k system with Gromacs 4.0.5)
# - The job reads from and writes in the directory in which it was launched. Used 
#   'initialdir = STARTDIR' to change  this behaviour.


# Set this to the same value as walltime; it ensures that Gromacs shuts down in time
# to write continuation files.
WALL_HOURS=0.33


# suggested environment settings (only change after benchmarking)
export MP_EAGER_LIMIT=65536
export MP_SHARED_MEMORY=yes
export MEMORY_AFFINITY=MCM
export MP_TASK_AFFINITY=MCM


MPIRUN=poe
MDRUN=/usr/local/packages/gromacs/g_4.0.5/bin/mdrun_mpi


# run gromacs
$MPIRUN $MDRUN -v -deffnm md -maxh ${WALL_HOURS}  -cpi


#########
# Notes
#########
# For finding best PME/pp-node partitioning (comment out 'run gromacs')
# (g_tune_pme uses MPIRUN and MDRUN environment variables!)
# export MPIRUN MDRUN
# GTUNEPME=/usr/local/packages/gromacs/g_4.0.5/bin/g_tune_pme

# $GTUNEPME -np $NPROC -v -deffnm md


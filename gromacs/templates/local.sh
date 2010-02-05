#!/bin/sh
#
# Run a gromacs MD job on the local machine.
# This script is very simplistic and assumes that you can simply run
# gromacs with mpiexec (see below)

# set this to the same value as walltime; mdrun will stop cleanly
# at 0.99 * WALL_HOURS 
WALL_HOURS=0.33

NCORES=8

# deffnm line is possibly modified by gromacs.setup
# (leave it as it is in the template)
DEFFNM=md

TPR=${DEFFNM}.tpr
OUTPUT=${DEFFNM}.out
PDB=${DEFFNM}.pdb

MDRUN_OPTS=""

mpiexec -n $NCORES mdrun_mpi -nice 19 -v -deffnm ${DEFFNM} -c ${PDB} -cpi \
    $MDRUN_OPTS \
    -maxh ${WALL_HOURS} > $OUTPUT

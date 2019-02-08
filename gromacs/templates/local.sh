#!/bin/sh
#
# Run a gromacs MD job on the local machine -- THIS IS JUST AN EXAMPLE
# This script is very simplistic and assumes that you can simply run
# gromacs with mpiexec (see below)

NCORES=8

# deffnm line is possibly modified by gromacs.setup
# (leave it as it is in the template)
DEFFNM=md

TPR=${DEFFNM}.tpr
OUTPUT=${DEFFNM}.out
PDB=${DEFFNM}.pdb

MDRUN_OPTS=""

mpiexec -n $NCORES mdrun_mpi -nice 19 -v -deffnm ${DEFFNM} -c ${PDB} -cpi -append \
    $MDRUN_OPTS >$OUTPUT 2>&1


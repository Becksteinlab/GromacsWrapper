#!/bin/sh
#
# Run a gromacs MD job on the local machine.
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

MPIRUN="mpiexec -n $NCORES"
MDRUN=mdrun_openmp64

$MPIRUN $MDRUN -nice 19 -v -deffnm ${DEFFNM} -c ${PDB} -cpi \
    $MDRUN_OPTS > $OUTPUT
rc=$?


exit $rc
#!/usr/bin/env bash 

#set -e

# Clean
rm \#* confout.gro *.log output.top mdout.mdp topol.tpr traj.trr ener.edr *.xvg *xtc *.cpt *~

# Reference run
grompp -f ../grompp.mdp -p processed.top -maxwarn 1
mdrun -v

# Rerun
mdrun -v -rerun traj.xtc
mv md.log desired.log
g_energy -o desired.xvg << EOF
Bond
U-B
Proper-Dih.
Improper-Dih.
CMAP-Dih.
LJ-14
LJ-(SR)
Coulomb-14
Coulomb-(SR)
Coul.-recip.

EOF

python ../test.py

grompp -f ../grompp.mdp -p output.top -maxwarn 1
mdrun -rerun traj.xtc -v 
mv md.log actual.log
g_energy -o actual.xvg << EOF
Bond
U-B
Proper-Dih.
Improper-Dih.
CMAP-Dih.
LJ-14
LJ-(SR)
Coulomb-14
Coulomb-(SR)
Coul.-recip.

EOF

python ../test_energies.py

rm \#*




from gromacs.fileformats import TOP
import numpy as np
import math
import copy, argparse


def scale_angles(mol, angles):
	new_angles = {}
	for dh in mol.angles:
		atypes = dh.atom1.get_atomtype(), dh.atom2.get_atomtype(), dh.atom3.get_atomtype()
		atypes = [a.replace("_", "").replace("=","") for a in atypes]
		for iswitch in range(16):
			if  (iswitch%2==0 ):
				a1=atypes[0]; a2=atypes[1]; a3=atypes[2]
			else:
				a1=atypes[2]; a2=atypes[1]; a3=atypes[0]

			if((iswitch//2)%2==1): a1="X";
			if((iswitch//4)%2==1): a2="X";
			if((iswitch//8)%2==1): a3="X";
			key = "{0}-{1}-{2}-{3}".format(a1, a2, a3, dh.gromacs['func'])
			if (key in angles): 
				for i, at in enumerate(angles[key]):
					#new_angles.append(at)
					new_angles[key] = at
				break
	return new_angles.values()


def scale_dihedrals(mol, dihedrals):
	new_dihedrals = {}
	for dh in mol.dihedrals:
		atypes = dh.atom1.get_atomtype(), dh.atom2.get_atomtype(), dh.atom3.get_atomtype(), dh.atom4.get_atomtype()
		atypes = [a.replace("_", "").replace("=","") for a in atypes]
		for iswitch in range(32):
			if  (iswitch%2==0 ):
				a1=atypes[0]; a2=atypes[1]; a3=atypes[2]; a4=atypes[3]
			else:
				a1=atypes[3]; a2=atypes[2]; a3=atypes[1]; a4=atypes[0]

			if((iswitch//2)%2==1): a1="X";
			if((iswitch//4)%2==1): a2="X";
			if((iswitch//8)%2==1): a3="X";
			if((iswitch//16)%2==1): a4="X";
			key = "{0}-{1}-{2}-{3}-{4}".format(a1, a2, a3, a4, dh.gromacs['func'])
			if (key in dihedrals): 
				for i, dt in enumerate(dihedrals[key]):
					#new_dihedrals.append(dt)
					new_dihedrals[key] = dt
				break

	print new_dihedrals
	return new_dihedrals.values()


def scale_impropers(mol, impropers):
	new_impropers = {}
	for im in mol.impropers:
		atypes = im.atom1.get_atomtype(), im.atom2.get_atomtype(), im.atom3.get_atomtype(), im.atom4.get_atomtype()
		atypes = [a.replace("_", "").replace("=","") for a in atypes]
		for iswitch in range(32):
			if  (iswitch%2==0 ):
				a1=atypes[0]; a2=atypes[1]; a3=atypes[2]; a4=atypes[3];
			else:
				a1=atypes[3]; a2=atypes[2]; a3=atypes[1]; a4=atypes[0];
			if((iswitch/2)%2==1): a1="X";
			if((iswitch/4)%2==1): a2="X";
			if((iswitch/8)%2==1): a3="X";
			if((iswitch/16)%2==1): a4="X";
			key = "{0}-{1}-{2}-{3}-{4}".format(a1, a2, a3, a4, im.gromacs['func'])
			if (key in impropers): 
				for i, imt in enumerate(impropers[key]):
					new_impropers[key] = imt
				break
	print new_impropers
	return new_impropers.values()


parser = argparse.ArgumentParser()
parser.add_argument("input")
parser.add_argument("output")
args = parser.parse_args()

top = TOP(args.input)
molname = top.molecules[0].name
mol = top.dict_molname_mol[molname]

#
# ATOMTYPES
#
atomtypes = set([a.atomtype for a in mol.atoms])
top.atomtypes = [at for at in top.atomtypes if at.atype in atomtypes]

#
# BONDTYPES
#
bondtypes = set([ tuple(sorted((b.atom1.atomtype, b.atom2.atomtype))) for b in mol.bonds])
bondtypes_dictionary = {tuple(sorted((bt.atype1, bt.atype2))): bt for bt in top.bondtypes}
top.bondtypes = [bondtypes_dictionary[bt] for bt in bondtypes]

#
# Build bond dictionary
# 
angletypes = {}
for at in top.angletypes:
	name = "{0}-{1}-{2}-{3}".format(at.atype1, at.atype2, at.atype3, at.gromacs['func'])
	if not name in angletypes: angletypes[name] = []
	angletypes[name].append(at)

#
# Build dihedral dictionary
#
dihedraltypes = {}
for dt in top.dihedraltypes:
	name = "{0}-{1}-{2}-{3}-{4}".format(dt.atype1, dt.atype2, dt.atype3, dt.atype4, dt.gromacs['func'])
	if not name in dihedraltypes: dihedraltypes[name] = []
	dihedraltypes[name].append(dt)
print("Build dihedraltypes dictionary with {0} entries".format(len(dihedraltypes)))

#
# Build improper dictionary
#
impropertypes = {}
for it in top.impropertypes:
	name = "{0}-{1}-{2}-{3}-{4}".format(it.atype1, it.atype2, it.atype3, it.atype4, it.gromacs['func'])
	if not name in impropertypes: impropertypes[name] = []
	impropertypes[name].append(it)
print("Build impropertypes dictionary with {0} entries".format(len(impropertypes)))

top.angletypes = scale_angles(mol, angletypes)
top.dihedraltypes = scale_dihedrals(mol, dihedraltypes)
top.impropertypes = scale_impropers(mol, impropertypes)

top.nonbond_params = []
top.cmaptypes = []
atomtypes = set([at.atype for at in top.atomtypes])

pairtypes = [pt for pt in top.pairtypes if (pt.atype1 in atomtypes) and (pt.atype2 in atomtypes)]
top.pairtypes = pairtypes

# Remove non-default moleculetypes
for k in top.dict_molname_mol.keys():
	if k in [molname]: continue
	del top.dict_molname_mol[k]

top.write(args.output)

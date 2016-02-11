

#import pytopol.parsers.grotop
#from pytopol.parsers.grotop import SystemToGroTop

from gromacs.fileformats import TOP

def test_defaults(actual, desired):
	print("defaults (not implemented)")

def test_atomtypes(actual, desired):
	at1 = [(at.atype, at.mass, at.charge, at.gromacs) for at in actual.atomtypes]
	at2 = [(at.atype, at.mass, at.charge, at.gromacs) for at in desired.atomtypes]
	assert len(at1) == len(at2), "number of atomtypes is not identical"
	assert at1 == at2, "terms in atomtypes are not identical"
	print("atomtypes ok")

def test_pairtypes(actual, desired):
	pt1 = [(pt.atype1, pt.atype2, pt.gromacs) for pt in actual.pairtypes]
	pt2 = [(pt.atype1, pt.atype2, pt.gromacs) for pt in desired.pairtypes]
	assert len(pt1) == len(pt2), "number of pairtypes is not identical"
	assert pt1 == pt2, "terms in pairtypes are not identical"
	print("pairtypes ok")

def test_bondtypes(actual, desired):
	bt1 = [(bt.atype1, bt.atype2, bt.gromacs) for bt in actual.pairtypes]
	bt2 = [(bt.atype1, bt.atype2, bt.gromacs) for bt in desired.pairtypes]
	assert len(bt1) == len(bt2), "number of bondtypes is not identical"
	assert bt1 == bt2, "terms in bondtypes are not identical"
	print("bondtypes ok")

def test_constrainttypes(actual, desired):
	ct1 = [(ct.atype1, ct.atype2, ct.gromacs) for ct in actual.constrainttypes]
	ct2 = [(ct.atype1, ct.atype2, ct.gromacs) for ct in desired.constrainttypes]
	assert len(ct1) == len(ct2), "number of constrainttypes is not identical"
	assert ct1 == ct2, "terms in constrainttypes are not identical"
	print("constrainttypes ok")

def test_angletypes(actual, desired):
	at1 = [(at.atype1, at.atype2, at.atype3, at.gromacs) for at in actual.angletypes]
	at2 = [(at.atype1, at.atype2, at.atype3, at.gromacs) for at in desired.angletypes]
	assert len(at1) == len(at2), "number of angletypes is not identical"
	assert at1 == at2, "terms in angletypes are not identical"
	print("angletypes ok")
	
def test_dihedraltypes(actual, desired):
	dt1 = [(dt.atype1, dt.atype2, dt.atype3, dt.atype4, dt.gromacs) for dt in actual.dihedraltypes]
	dt2 = [(dt.atype1, dt.atype2, dt.atype3, dt.atype4, dt.gromacs) for dt in desired.dihedraltypes]
	assert len(dt1) == len(dt2), "number of dihedraltypes is not identical"
	assert dt1 == dt2, "terms in dihedraltypes are not identical"
	print("dihedraltypes ok")

def test_cmaptypes(actual, desired):
	ct1 = [(ct.atype1, ct.atype2, ct.atype3, ct.atype4, ct.atype5, ct.atype6, ct.atype7, ct.atype8, ct.gromacs) for ct in actual.cmaptypes]
	ct2 = [(ct.atype1, ct.atype2, ct.atype3, ct.atype4, ct.atype5, ct.atype6, ct.atype7, ct.atype8, ct.gromacs) for ct in desired.cmaptypes]
	assert len(ct1) == len(ct2), "number of cmaptypes is not identical"
	assert ct1 == ct2, "terms in cmaptypes are not identical"
	print("cmaptypes ok")

def test_moleculetypes(actual, desired):
	print("moleculetypes (not implemented)")

def test_molecules(actual, desired):
	mol1 = [(mol.name, len(mol.atoms), len(mol.bonds), len(mol.angles), len(mol.dihedrals)) for mol in  actual.molecules]
	mol2 = [(mol.name, len(mol.atoms), len(mol.bonds), len(mol.angles), len(mol.dihedrals)) for mol in  desired.molecules]
	assert len(mol1) == len(mol2), "number of molecules is not identical"
	assert mol1 == mol2, "number of terms in molecules are not identical"
	print("molecules ok")


def test_top(actual, desired):
	# TODO test_defaults
	test_defaults(actual, desired)
	test_atomtypes(actual, desired)
	test_pairtypes(actual, desired)
	test_constrainttypes(actual, desired)

	# Angles, dihedrals and CMAP
	test_angletypes(actual, desired)
	test_dihedraltypes(actual, desired)
	test_cmaptypes(actual, desired)

	# TODO test moleculetypes	
	test_moleculetypes(actual, desired)
	test_molecules(actual, desired)

desired = TOP("processed.top")
desired.write("output.top")
actual = TOP("output.top")


# Weak test but if that's not true, there is no point in contining
assert str(actual) == str(desired)

test_top(actual, desired)
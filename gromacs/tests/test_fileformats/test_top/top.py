from __future__ import division, absolute_import, print_function

import os
import numpy as np
import numpy.matlib
from numpy.testing import assert_array_equal, assert_, run_module_suite

from gromacs.fileformats import TOP
from gromacs import testing as tm

def helper(attr, attr1, attr2):
	print(attr)
	for pair in zip(attr1, attr2): 
		if pair[0] == pair[1]: continue
		print(pair)

class TopologyTest(object):

	def test_basic(self):
	    path = tm.get_data_path(self.processed)
	    top = TOP(path)
	    assert_(top.dict_molname_mol.keys() == self.molecules)


	def test_equal(self):
		"""
			Load the same topology twice and check if __eq__ comparisions work
		"""
		path = tm.get_data_path(self.processed)

		top1 = TOP(path)
		top2 = TOP(path)

		attrs1 = [section for section in top1.found_sections if "types" in section]
		attrs2 = [section for section in top1.found_sections if "types" in section]

		assert_(attrs1==attrs2)

		for attr in attrs1:
			assert_(getattr(top1, attr) == getattr(top2, attr), attr)

	def test_parameter_types(self):
		"""Test if all the parameter types are the same across two topologies
		"""
		pass

	def test_parameters(self):
		"""Test if per-molecule parameters and parameter type assignments are identical
		"""
		pass

	def test_molecule_parameters(self):
		"""Called by `test_parameters()` for each molecule in the system
		"""

	def test_read_write(self):
		"""Read a topology, write it out, and read in the output again.
		Writing the topology out should make no change to the topology. 
		"""
		path = tm.get_data_path(self.processed)
		filename = '/tmp/processed-%s.top' % os.getpid()

		top1 = TOP(path)
		top1.write(filename)

		# make life harder, write out again
		top2 = TOP(filename)
		top2.write(filename)

		top2 = TOP(filename)

		attrs1 = [section for section in top1.found_sections if "types" in section]
		attrs2 = [section for section in top1.found_sections if "types" in section]

		assert_(attrs1==attrs2)
		#attrs1 = ['atomtypes', 'pairtypes', 'bondtypes', 'constrainttypes', 'angletypes', 'dihedraltypes', 'dihedraltypes', 'cmaptypes']
		#attrs1 = ['atomtypes',]

		for attr in attrs1:
			attr1 = getattr(top1, attr)
			attr2 = getattr(top2, attr)
			assert_(attr1 == attr2, helper(attr, attr1, attr2))


if __name__ == "__main__":
    run_module_suite()
from __future__ import division, absolute_import, print_function

import numpy as np
import numpy.matlib
from numpy.testing import assert_array_equal, assert_, run_module_suite

from gromacs.fileformats import TOP
from gromacs import testing as tm

class TopologyPrimitive(object):

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

		#attrs = "angletypes","cmaptypes","defaults","fname","impropertypes","nonbond_params","atomtypes","constrainttypes","pairtypes","bondtypes","dihedraltypes","interactiontypes"

		for attr in attrs1:
			assert_(getattr(top1, attr) == getattr(top2, attr), attr)


if __name__ == "__main__":
    run_module_suite()
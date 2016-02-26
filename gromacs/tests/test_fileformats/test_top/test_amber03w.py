
# GromacsWrapper: test_amber03w.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

from __future__ import division, absolute_import, print_function

import numpy as np
import numpy.matlib
from numpy.testing import assert_array_equal, assert_, run_module_suite

from gromacs.fileformats import TOP
from gromacs import testing as tm

from top import TopologyTest

class TestCharmm(TopologyTest):
	processed = 'amber03w/processed.top'
	conf = 'amber03w/conf.gro'
	molecules = ['Protein_chain_A', 'SOL', 'IB+', 'CA', 'CL', 'NA', 'MG', 'K', 'RB', 'CS', 'LI', 'ZN']

if __name__ == "__main__":
    run_module_suite()

# GromacsWrapper: test_basic.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

from __future__ import division, absolute_import, print_function

import numpy as np
import numpy.matlib
from numpy.testing import assert_array_equal, assert_, run_module_suite

def test_empty():
    assert_(2+4 == 6)



if __name__ == "__main__":
    run_module_suite()
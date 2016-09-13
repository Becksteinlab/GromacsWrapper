# GromacsWrapper: test_example.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

from __future__ import division, absolute_import, print_function

import pytest

from numpy.testing import assert_almost_equal

from gromacs import cbook
import gromacs.setup

from gromacs.tests.datafiles import datafile

def test_grompp_qtot(tmpdir):
    pdb = datafile("1ake_A.pdb")
    top = tmpdir.mkdir("top")
    with top.as_cwd():
        f = gromacs.setup.topology(struct=pdb, ff="oplsaa", water="tip4p")
        with open('none.mdp','w') as mdp:
            mdp.write('; empty mdp file\nrcoulomb = 1\nrvdw = 1\nrlist = 1\n')
        qtot = cbook.grompp_qtot(f="none.mdp", c=f['struct'], p=f['top'],
                                 stdout=False, maxwarn=10)
    assert_almost_equal(qtot, -4, decimal=5,
                        err_msg="grompp_qtot() failed to compute total charge correctly")

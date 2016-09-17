# GromacsWrapper: test_example.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

from __future__ import division, absolute_import, print_function

import os.path
import pytest

from numpy.testing import assert_almost_equal

import gromacs.setup

from gromacs.tests.datafiles import datafile

def test_trj_compact_main(tmpdir):
    pdb = datafile("1ake_A.pdb")
    top = tmpdir.mkdir("top")
    mdpfile = "simple.mdp"
    tprfile = "simple.tpr"
    outfile = "compact.pdb"
    with top.as_cwd():
        f = gromacs.setup.topology(struct=pdb, ff="oplsaa", water="tip4p")
        with open(mdpfile, 'w') as mdp:
            mdp.write('; empty mdp file\nrcoulomb = 1\nrvdw = 1\nrlist = 1\n')
        gromacs.grompp(f=mdpfile, o=tprfile, c=f["struct"], p=f["top"])
        gromacs.setup.trj_compact_main(s=tprfile, f=f["struct"], o=outfile, input=("protein", "system"))
        assert os.path.exists(outfile)

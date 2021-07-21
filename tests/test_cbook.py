# GromacsWrapper: test_example.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

from __future__ import division, absolute_import, print_function

import os.path
import re

import pytest
from numpy.testing import assert_almost_equal

from gromacs import cbook
import gromacs.setup

from .datafiles import datafile


@pytest.fixture
def simulation(tmpdir):
    pdb = datafile("1ake_A.pdb")
    with tmpdir.mkdir("top").as_cwd():
        f = gromacs.setup.topology(struct=pdb, ff="oplsaa", water="tip4p")
        yield f

def test_grompp_qtot(tmpdir, simulation):
    with tmpdir.mkdir("qtot").as_cwd():
        with open('none.mdp', 'w') as mdp:
            mdp.write('; empty mdp file\nrcoulomb = 1\nrvdw = 1\nrlist = 1\n')
        qtot = cbook.grompp_qtot(f="none.mdp", c=simulation['struct'], p=simulation['top'],
                                 stdout=False, maxwarn=10)
    assert_almost_equal(qtot, -4, decimal=5,
                        err_msg="grompp_qtot() failed to compute total charge correctly")


def test_portable_topology(tmpdir, simulation):
    with tmpdir.mkdir("processed").as_cwd():
        pptopol = cbook.create_portable_topology(simulation['top'], simulation['struct'])

    # correct filename
    assert os.path.split(pptopol)[-1].startswith("pp_")

    lines =  open(pptopol).readlines()
    subsections_only = [line.strip() for line in lines if line.startswith("[")]

    assert re.match(r";\s+File 'system.top' was generated", lines[1])
    assert re.match(r";\s+This is a standalone topology file", lines[6])
    for subsection in ("[ defaults ]", "[ atomtypes ]", "[ bondtypes ]",
                       "[ constrainttypes ]", "[ angletypes ]", "[ dihedraltypes ]",
                       "[ moleculetype ]", "[ atoms ]", "[ bonds ]", "[ pairs ]",
                       "[ angles ]", "[ dihedrals ]",
                       "[ system ]", "[ molecules ]"):
        assert subsection in subsections_only




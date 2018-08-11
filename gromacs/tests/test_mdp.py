# GromacsWrapper: test_example.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

import gromacs
import pytest
from numpy.testing import assert_equal

from .datafiles import datafile

@pytest.fixture(
    params=['original', 'written'],
)
def CUSTOM_EM_MDP(request, tmpdir):
    mdp = gromacs.fileformats.mdp.MDP(datafile('custom_em.mdp'))
    if request.param == 'original':
        return mdp
    elif request.param == 'written':
        # to check that written mdp has same data
        out = str(tmpdir.join('out.mdp'))

        mdp.write(out)

        return gromacs.fileformats.mdp.MDP(out)

class TestMDP(object):
    def test_values(self, CUSTOM_EM_MDP):
        mdp = CUSTOM_EM_MDP
        assert_equal(mdp['include'], ['-I.', '-I..', '-I../top'])
        assert_equal(mdp['define'], ['-DPOSRES', '-DFLEXSPC'])
        assert mdp['integrator'] == 'cg'
        assert mdp['emtol'] == 500
        assert mdp['emstep'] == 0.01
        assert mdp['nsteps'] == 1000
        assert mdp['nstcgsteep'] == 100
        assert mdp['constraints'] == 'none'
        assert mdp['nstcomm'] == 1
        assert mdp['cutoff-scheme'] == 'Verlet'
        assert mdp['vdwtype'] == 'cutoff'
        assert mdp['coulombtype'] == 'PME'
        assert mdp['ns_type'] == 'grid'
        assert mdp['rlist'] == 1.0
        assert mdp['rcoulomb'] == 1.0
        assert mdp['rvdw'] == 1.0
        assert mdp['rvdw-switch'] == 0.8
        assert mdp['Tcoupl'] == 'no'
        assert mdp['Pcoupl'] == 'no'
        assert mdp['gen_vel'] == 'no'
        assert mdp['nstxout'] == 0

    def test_comments(self, CUSTOM_EM_MDP):
        mdp = CUSTOM_EM_MDP

        assert mdp['C0001'] == 'custom EM'
        assert sum(1 for k in mdp if k.startswith('C0')) == 1

    def test_blank_lines(self, CUSTOM_EM_MDP):
        mdp = CUSTOM_EM_MDP
        assert sum(1 for k in mdp if k.startswith('B0')) == 6

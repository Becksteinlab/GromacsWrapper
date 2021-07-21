# GromacsWrapper: test_example.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

import gromacs
import pytest
from numpy.testing import assert_equal

from ..datafiles import datafile

@pytest.fixture(
    params=['original', 'written', 'no_autoconvert'],
)
def CUSTOM_EM_MDP(request, tmpdir):
    autoconvert = not request.param == 'no_autoconvert'
    mdp = gromacs.fileformats.mdp.MDP(datafile('custom_em.mdp'),
                                      autoconvert=autoconvert)
    if request.param == 'written':
        # to check that written mdp has same data
        out = str(tmpdir.join('out.mdp'))
        mdp.write(out)
        mdp = gromacs.fileformats.mdp.MDP(out)
    return mdp, autoconvert


class TestMDP(object):
    def test_values(self, CUSTOM_EM_MDP):
        mdp, autoconvert = CUSTOM_EM_MDP

        if autoconvert:
            def conv(val):
                return val
        else:
            def conv(val):
                return str(val)

        if autoconvert:
            assert_equal(mdp['include'], ['-I.', '-I..', '-I../top'])
            assert_equal(mdp['define'], '-DPOSRES')
        else:
            assert mdp['include'] == '-I. -I.. -I../top'
            assert mdp['define'] == '-DPOSRES'
        assert mdp['integrator'] == 'cg'
        assert mdp['emtol'] == conv(500)
        assert mdp['emstep'] == conv(0.01)
        assert mdp['nsteps'] == conv(1000)
        assert mdp['nstcgsteep'] == conv(100)
        assert mdp['constraints'] == 'none'
        assert mdp['nstcomm'] == conv(1)
        assert mdp['cutoff-scheme'] == 'Verlet'
        assert mdp['vdwtype'] == 'cutoff'
        assert mdp['coulombtype'] == 'PME'
        assert mdp['ns_type'] == 'grid'
        assert mdp['rlist'] == conv(1.0)
        assert mdp['rcoulomb'] == conv(1.0)
        assert mdp['rvdw'] == conv(1.0)
        assert mdp['rvdw-switch'] == conv(0.8)
        assert mdp['Tcoupl'] == 'no'
        assert mdp['Pcoupl'] == 'no'
        assert mdp['gen_vel'] == 'no'
        assert mdp['nstxout'] == conv(0)

    def test_comments(self, CUSTOM_EM_MDP):
        mdp, _ = CUSTOM_EM_MDP

        assert mdp['C0001'] == 'custom EM'
        assert sum(1 for k in mdp if k.startswith('C0')) == 1

    def test_blank_lines(self, CUSTOM_EM_MDP):
        mdp, _ = CUSTOM_EM_MDP
        assert sum(1 for k in mdp if k.startswith('B0')) == 6

    def test_no_filename(self, tmpdir):
        mdp = gromacs.fileformats.mdp.MDP()

        mdp['define'] = ['-DPOSRES', '-DTHIS']
        mdp['vdwtype'] = 'cutoff'
        mdp['nsteps'] = 1234

        out = str(tmpdir.join('out.mdp'))
        mdp.write(out)

        back = gromacs.fileformats.mdp.MDP(out)

        assert_equal(back['define'], ['-DPOSRES', '-DTHIS'])
        assert back['vdwtype'] == 'cutoff'
        assert back['nsteps'] == 1234

@pytest.fixture
def NONSENSE_MDP(tmpdir):
    outfile = str(tmpdir.join('nonsense.mdp'))

    with open(datafile('custom_em.mdp'), 'r') as infile:
        data = infile.read()

    data += 'errors: plenty\n'

    with open(outfile, 'w') as out:
        out.write(data)

    return outfile


def test_bad_mdp(NONSENSE_MDP):
    with pytest.raises(gromacs.ParseError,
                       match="unknown line in mdp file, 'errors: plenty'"):
        gromacs.fileformats.mdp.MDP(NONSENSE_MDP)

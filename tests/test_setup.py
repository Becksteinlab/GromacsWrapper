# GromacsWrapper: test_example.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

from __future__ import division, absolute_import, print_function

import os.path
import pytest

from numpy.testing import assert_almost_equal

import gromacs.setup

from .datafiles import datafile

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
        gromacs.setup.trj_compact_main(s=tprfile, f=f["struct"], o=outfile,
                                       input=("protein", "system"))
        assert os.path.exists(outfile)

@pytest.fixture(scope="session")
def topology(tmpdir_factory, struct=datafile("1ake_A_protein.pdb")):
    # note: use protein-only input 1ake_A_protein.pdb because solvation fails
    #       if crystal waters are included (in 1ake_A.pdb)
    TMPDIR = tmpdir_factory.mktemp('1ake')
    with TMPDIR.as_cwd():
        topol_args = gromacs.setup.topology(struct=struct, ff="oplsaa", water="tip4p")
    return TMPDIR, topol_args

@pytest.fixture(scope="session")
def solvate(topology):
    TMPDIR, topol_args = topology
    with TMPDIR.as_cwd():
        solvate_args = gromacs.setup.solvate(concentration=0.15,
                                             water="tip4p",
                                             **topol_args)
    return TMPDIR, solvate_args

@pytest.fixture(scope="session")
def energy_minimize(solvate, low_performance):
    # run short energy minimization with cheapest minimizer
    TMPDIR, solvate_args = solvate
    nt = 2 if low_performance else 0
    with TMPDIR.as_cwd():
        em_args = gromacs.setup.energy_minimize(mdrun_args={'nt': nt},
                                                integrator="steep",
                                                emtol=5000,
                                                maxwarn=1,
                                                **solvate_args)
    return TMPDIR, em_args


def test_topology(topology):
    TMPDIR, topol_args = topology
    top = topol_args['top']
    struct = topol_args['struct']
    posres = topol_args['posres']
    assert os.path.exists(top)
    assert os.path.exists(struct)
    assert os.path.exists(posres)
    # add more tests for content of files!

def test_solvate(solvate):
    TMPDIR, solvate_args = solvate

    assert_almost_equal(solvate_args['qtot'], 0.0)
    assert os.path.exists(solvate_args['struct'])
    assert os.path.exists(solvate_args['ndx'])
    assert solvate_args['mainselection'] == '"Protein"'
    # add more tests for content of files!

def test_energy_minimize(energy_minimize):
    TMPDIR, em_args = energy_minimize
    assert os.path.exists(em_args['struct'])
    assert os.path.exists(em_args['top'])
    assert em_args['mainselection'] == '"Protein"'
    # add more tests for content of files!

def test_energy_minimize_custom_mdp(solvate, low_performance,
                                    mdp=datafile("custom_em.mdp")):
    TMPDIR, solvate_args = solvate
    nt = 2 if low_performance else 0
    with TMPDIR.as_cwd():
        try:
            em_args = gromacs.setup.energy_minimize(mdrun_args={'nt': nt},
                                                    mdp=mdp,
                                                    emtol=5000,
                                                    **solvate_args)
        except gromacs.exceptions.GromacsError as err:
            # sometimes the em does not converge at all, e.g. 5.02988e+04 on atom 3277;
            # (happens on Travis Linux with Gromacs 4.6.5 but not locally or on Travis OSX) so we
            # re-run with a ridiculous tolerance so that we can at least test that the whole
            # function can run to completion
            em_args = gromacs.setup.energy_minimize(mdrun_args={'nt': nt},
                                                    mdp=mdp,
                                                    emtol=6e4,
                                                    **solvate_args)
    assert os.path.exists(em_args['struct'])
    assert os.path.exists(em_args['top'])
    assert em_args['mainselection'] == '"Protein"'
    # add more tests for content of files!

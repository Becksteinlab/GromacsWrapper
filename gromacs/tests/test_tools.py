# GromacsWrapper: test_example.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

from __future__ import division, absolute_import, print_function

import pytest

import gromacs

common_tool_names = ["pdb2gmx", "grompp", "editconf", "mdrun"]
aliased_tool_names = list(gromacs.tools.NAMES5TO4.values())

@pytest.fixture(scope="function",
                params=set(common_tool_names + aliased_tool_names))
def gromacs_tool(request):
    return getattr(gromacs, request.param)

def test_tools_help(gromacs_tool):
    rc, out, err = gromacs_tool(h=True, stdout=False, stderr=False)
    assert rc == 0, "Gromacs command {0} failed".format(gromacs_tool.command_name)
    assert out + err, "Gromacs command {0} produced no output for -h".format(
        gromacs_tool.command_name)

def test_failure_raises():
    # unknown option
    with pytest.raises(gromacs.GromacsError):
        gromacs.grompp(y=True)

def test_failure_warns():
    # unknown option
    grompp_warn = gromacs.tools.Grompp(failure="warn")
    with pytest.warns(gromacs.GromacsFailureWarning):
        grompp_warn(y=True)

def test_failure_ignore():
    # unknown option
    grompp_ignore = gromacs.tools.Grompp(failure=None)
    try:
        grompp_ignore(y=True)
    except Exception as err:
        raise AssertionError("Should have ignored exception {}".format(err))

import pathlib

path_gmx = pathlib.Path('/usr/local/bin/gmx')
path_gmx_mpi = pathlib.Path('~/gmx_mpi').expanduser()
path_gmx_mpi.symlink_to(str(path_gmx))

gromacs.config.cfg['Gromacs']['tools'] = str(path_gmx_mpi)
with open('/Users/theavey/.gromacswrapper.cfg', 'w') as cfg_fo:
    gromacs.config.cfg.write(cfg_fo)

import gromacs as new_gromacs



@pytest.fixture(scope="function",
                params=set(common_tool_names + aliased_tool_names))
def new_gromacs_tool(request):
    return getattr(new_gromacs, request.param)

def test_tools_help(gromacs_tool):
    rc, out, err = gromacs_tool(h=True, stdout=False, stderr=False)
    assert rc == 0, "Gromacs command {0} failed".format(gromacs_tool.command_name)
    assert out + err, "Gromacs command {0} produced no output for -h".format(
        gromacs_tool.command_name)

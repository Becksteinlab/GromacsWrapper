# GromacsWrapper: test_example.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

from __future__ import division, absolute_import, print_function

import pytest

import gromacs

common_tool_names = ["pdb2gmx", "grompp", "editconf", "mdrun"]
aliased_tool_names = gromacs.tools.NAMES5TO4.values()

@pytest.fixture(scope="module",
                params=set(common_tool_names + aliased_tool_names))
def gromacs_tool(request):
    return getattr(gromacs, request.param)

def test_tools_help(gromacs_tool):
    rc, out, err = gromacs_tool(h=True, stdout=False, stderr=False)
    assert rc == 0, "Gromacs command {0} failed".format(gromacs_tool.command_name)
    assert out + err, "Gromacs command {0} produced no output for -h".format(
        gromacs_tool.command_name)



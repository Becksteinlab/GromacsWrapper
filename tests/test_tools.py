# GromacsWrapper: test_example.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

import pytest

import gromacs


common_tool_names = ["pdb2gmx", "grompp", "editconf", "mdrun"]
aliased_tool_names = list(gromacs.tools.NAMES5TO4.values())

changed_tools = {
    "do_dssp": "dssp",  # GMX 2022.x --> 2023.x
}


@pytest.fixture(scope="function", params=set(common_tool_names + aliased_tool_names))
def gromacs_tool(request):
    # for testing, work around changed tools that cannot be aliased
    # via NAMES5TO4 because the call signature has changed: replace
    # old toolname with the functionally equivalent new name
    try:
        tool = getattr(gromacs, request.param)
    except AttributeError:
        tool = getattr(gromacs, changed_tools[request.param])
    return tool


def test_tools_help(gromacs_tool):
    rc, out, err = gromacs_tool(h=True, stdout=False, stderr=False)
    assert rc == 0, "Gromacs command {0} failed".format(gromacs_tool.command_name)
    assert out + err, "Gromacs command {0} produced no output for -h".format(
        gromacs_tool.command_name
    )


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


class TestRelease(object):
    # add tested releases here
    major_releases = (
        "4",
        "5",
        "2016",
        "2018",
        "2019",
        "2020",
        "2021",
        "2022",
        "2023",
        "2024",
    )

    def test_release(self):
        assert gromacs.release().startswith(self.major_releases)

    def test_release_startswith(self):
        assert gromacs.release.startswith(self.major_releases)

    def test_str(self):
        assert str(gromacs.release()) == gromacs.release()

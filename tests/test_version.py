# GromacsWrapper: test_example.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

import re
import pytest

import gromacs


def test_version():
    release = gromacs.__version__
    assert isinstance(release, str)


def test_version_pep440_like():
    release = gromacs.__version__
    if release == "1+unknown":
        # testing in CI where versioningit does not have enough information to
        # generate the proper version information
        return
    match = re.match(r"\d+\.\d+\.\d+", release)
    assert match, f"Version {release} does not look like MAJOR.MINOR.PATCH..."

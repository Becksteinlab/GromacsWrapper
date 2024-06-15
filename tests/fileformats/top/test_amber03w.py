# GromacsWrapper: test_amber03w.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

from __future__ import division, absolute_import, print_function

import pytest

import gromacs

from .top import TopologyTest
from ...datafiles import datafile


class TestAmber03w(TopologyTest):
    processed = datafile("fileformats/top/amber03w/processed.top")
    conf = datafile("fileformats/top/amber03w/conf.gro")
    molecules = [
        "Protein_chain_A",
        "SOL",
        "IB+",
        "CA",
        "CL",
        "NA",
        "MG",
        "K",
        "RB",
        "CS",
        "LI",
        "ZN",
    ]

    @pytest.mark.xfail(
        gromacs.release().startswith(("2022", "2023", "2024")),
        reason="issue #236 https://github.com/Becksteinlab/GromacsWrapper/issues/236",
    )
    def test_mdrun(self, tmpdir, low_performance):
        super(TestAmber03w, self).test_mdrun(tmpdir, low_performance)

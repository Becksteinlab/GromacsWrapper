# GromacsWrapper: test_amber03star.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

from __future__ import division, absolute_import, print_function

import pytest

import gromacs

from .top import TopologyTest
from ...datafiles import datafile

@pytest.mark.xfail(gromacs.release().startswith("2022"),
                   reason="issue https://github.com/Becksteinlab/GromacsWrapper/issues/236")
class TestAmber03star(TopologyTest):
        processed = datafile('fileformats/top/amber03star/processed.top')
        conf = datafile('fileformats/top/amber03star/conf.gro')
        molecules = ['Protein', 'SOL', 'IB+', 'CA', 'CL', 'NA', 'MG', 'K', 'RB', 'CS', 'LI', 'ZN']

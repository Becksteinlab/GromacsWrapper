# GromacsWrapper: test_amber03star.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.



import pytest

from gromacs.exceptions import GromacsError

from .top import TopologyTest
from ...datafiles import datafile

class TestAmber03star(TopologyTest):
        processed = datafile('fileformats/top/amber03star/processed.top')
        conf = datafile('fileformats/top/amber03star/conf.gro')
        molecules = ['Protein', 'SOL', 'IB+', 'CA', 'CL', 'NA', 'MG', 'K', 'RB', 'CS', 'LI', 'ZN']

        @pytest.mark.xfail(raises=ValueError, reason="Not currently maintained. See #61.")
        def test_read_write(self, tmpdir):
                super(TestAmber03star, self).test_read_write(tmpdir)

        @pytest.mark.xfail(raises=GromacsError, reason="Not currently maintained. See #61.")
        def test_mdrun(self, tmpdir):
                super(TestAmber03star, self).test_mdrun(tmpdir)

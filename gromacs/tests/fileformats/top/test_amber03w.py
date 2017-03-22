# GromacsWrapper: test_amber03w.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.



import pytest

from gromacs.exceptions import GromacsError

from .top import TopologyTest
from ...datafiles import datafile

class TestAmber03w(TopologyTest):
        processed = datafile('fileformats/top/amber03w/processed.top')
        conf = datafile('fileformats/top/amber03w/conf.gro')
        molecules = ['Protein_chain_A', 'SOL', 'IB+', 'CA', 'CL', 'NA', 'MG', 'K', 'RB', 'CS', 'LI', 'ZN']

        @pytest.mark.xfail(raises=ValueError, reason="Not currently maintained. See #61.")
        def test_read_write(self, tmpdir):
                super(TestAmber03w, self).test_read_write(tmpdir)

        @pytest.mark.xfail(raises=GromacsError, reason="Not currently maintained. See #61.")
        def test_grompp(self, tmpdir):
                super(TestAmber03w, self).test_grompp(tmpdir)

        @pytest.mark.xfail(raises=GromacsError, reason="Not currently maintained. See #61.")
        def test_mdrun(self, tmpdir):
                super(TestAmber03w, self).test_mdrun(tmpdir)




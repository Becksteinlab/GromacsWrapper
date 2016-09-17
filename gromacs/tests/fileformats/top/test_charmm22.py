# GromacsWrapper: test_charmm22st.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

from __future__ import division, absolute_import, print_function

import pytest

from .top import TopologyTest
from ...datafiles import datafile

class TestCharmm22st(TopologyTest):
        processed = datafile('fileformats/top/charmm22st/processed.top')
        conf = datafile('fileformats/top/charmm22st/conf.gro')
        molecules = ['SOL', 'Protein', 'Ion', 'Cal', 'Ces', 'CL', 'K', 'NA', 'ZN']

        @pytest.mark.xfail(raises=ValueError, reason="Not currently maintained. See #61.")
        def test_read_write(self, tmpdir):
                super(TestCharmm22st, self).test_read_write(tmpdir)


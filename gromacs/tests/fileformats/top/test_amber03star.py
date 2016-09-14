# GromacsWrapper: test_amber03star.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

from __future__ import division, absolute_import, print_function

from .top import TopologyTest

class TestCharmm(TopologyTest):
        processed = 'fileformats/top/amber03star/processed.top'
        conf = 'fileformats/top/amber03star/conf.gro'
        molecules = ['Protein', 'SOL', 'IB+', 'CA', 'CL', 'NA', 'MG', 'K', 'RB', 'CS', 'LI', 'ZN']


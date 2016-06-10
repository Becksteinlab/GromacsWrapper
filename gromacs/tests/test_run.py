# GromacsWrapper: test_example.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

from __future__ import division, absolute_import, print_function

from .datafiles import datafile

import gromacs.run

class Test_check_mdrun_success(object):
    @staticmethod
    def test_no_logfile():
        assert(gromacs.run.check_mdrun_success("bogus_file.log") is None)

    @staticmethod
    def test_success_Gromacs4():
        assert(gromacs.run.check_mdrun_success(datafile('gromacs4_success.log')) is True)

    @staticmethod
    def test_incomplete_Gromacs4():
        assert(gromacs.run.check_mdrun_success(datafile('gromacs4_incomplete.log')) is False)

    @staticmethod
    def test_success_Gromacs5():
        assert(gromacs.run.check_mdrun_success(datafile('gromacs5_success.log')) is True)

    @staticmethod
    def test_incomplete_Gromacs5():
        assert(gromacs.run.check_mdrun_success(datafile('gromacs5_incomplete.log')) is False)




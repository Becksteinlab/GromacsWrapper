# GromacsWrapper: top.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

from __future__ import division, absolute_import, print_function

import os.path

import numpy as np

from numpy.testing import assert_array_equal
from pandas.util.testing import assert_frame_equal

import pytest

import gromacs
from gromacs.fileformats import TOP, XVG
from gromacs import scaling

from ...datafiles import datafile

def errmsg_helper(attr, attr1, attr2):
        msg = ["attr = " + str(attr)]
        for pair in zip(attr1, attr2):
                if pair[0] == pair[1]:
                        continue
                msg.append(str(pair))
        return "\n".join(msg)

def grompp(f, c, p, prefix="topol", **kwargs):
        s = prefix + '.tpr'
        po = prefix + '.mdp'

        rc, output, junk = gromacs.grompp(f=f, p=p, c=c, o=s, po=po, stdout=False, stderr=False,
                                          **kwargs)
        assert rc == 0, \
                "grompp -o {0} -po {1} -f {2} -c {3} -p {4} failed".format(s, po, f, c, p)
        return s

def mdrun(s, prefix, nt=0):
        o = prefix + '.trr'
        rc, output, junk = gromacs.mdrun(v=True, s=s, o=o, stdout=False, stderr=False, nt=nt)
        assert rc == 0, "mdrun failed"
        return o

def rerun_energy(s, o, prefix, nt=0):
        e = prefix + '.edr'
        rc, output, junk = gromacs.mdrun(v=True, s=s, rerun=o, e=e, stdout=False, stderr=False, nt=nt)
        assert rc == 0, "mdrun failed"

        xvg = prefix + '.xvg'
        rc, output, junk = gromacs.g_energy(f=e, o=xvg,
                                            input=("Proper-Dih.", "Improper-Dih.", "CMAP-Dih.",
                                                   "LJ-14", "Coulomb-14", "LJ-(SR)", "Coulomb-(SR)",
                                                   "Coul.-recip.", "Potential"), stdout=False, stderr=False)
        assert rc == 0, "g_energy failed"

        return XVG(xvg).to_df()

class TopologyTest(object):
        mdp = datafile('fileformats/top/grompp.mdp')

        def test_basic(self):
                path = self.processed
                top = TOP(path)
                assert list(top.dict_molname_mol.keys()) == self.molecules

        def test_equal(self):
                """
                        Load the same topology twice and check if __eq__ comparisions work
                """
                path = self.processed

                top1 = TOP(path)
                top2 = TOP(path)

                attrs1 = [section for section in top1.found_sections if "types" in section]
                attrs2 = [section for section in top1.found_sections if "types" in section]

                assert attrs1 == attrs2, "sections are not equal"

                for attr in attrs1:
                        assert getattr(top1, attr) == getattr(top2, attr), \
                                "{0} not identical".format(attr)

        # def test_parameter_types(self):
        #         """Test if all the parameter types are the same across two topologies
        #         """
        #         pass

        # def test_parameters(self):
        #         """Test if per-molecule parameters and parameter type assignments are identical
        #         """
        #         pass

        # def test_molecule_parameters(self):
        #         """Called by `test_parameters()` for each molecule in the system
        #         """

        def test_read(self, tmpdir):
                path = self.processed
                try:
                        top1 = TOP(path)
                except Exception as err:
                        raise AssertionError("Failed to read {0}. Raised:\n{1}".format(
                                path, str(err)))
                attrs1 = [section for section in top1.found_sections if "types" in section]
                assert attrs1

        def test_read_write(self, tmpdir):
                """Read a topology, write it out, and read in the output again.
                Writing the topology out should make no change to the topology.
                """
                path = self.processed
                with tmpdir.as_cwd():
                        filename1 = 'processed_1.top'
                        filename2 = 'processed_2.top'

                        top1 = TOP(path)
                        top1.write(filename1)

                        # make life harder, write out again
                        top2 = TOP(filename1)
                        top2.write(filename2)

                        top2 = TOP(filename2)

                        attrs1 = [section for section in top1.found_sections if "types" in section]
                        attrs2 = [section for section in top1.found_sections if "types" in section]

                        assert attrs1 == attrs2

                        #attrs1 = ['atomtypes', 'pairtypes', 'bondtypes', 'constrainttypes', 'angletypes', 'dihedraltypes', 'dihedraltypes', 'cmaptypes']
                        #attrs1 = ['atomtypes',]

                        for attr in attrs1:
                                attr1 = getattr(top1, attr)
                                attr2 = getattr(top2, attr)
                                assert attr1 == attr2, errmsg_helper(attr, attr1, attr2)

        def test_grompp(self, tmpdir):
                """Check if grompp can be run successfully at all"""
                f = self.mdp
                c = self.conf
                p = self.processed
                with tmpdir.as_cwd():
                        o = 'topol.tpr'
                        po = 'mdout.mdp'
                        rc, output, junk = gromacs.grompp(f=f, p=p, c=c, o=o, po=po, stdout=False, stderr=False)
                        assert rc == 0, "grompp -f {0} -o {1} ... failed to run".format(f, o)

        def test_mdrun(self, tmpdir, low_performance):
                """Check if grompp can be run successfully at all"""
                # set low_performance with
                #
                #    pytest --low-performance gromacs/tests
                #
                f = self.mdp
                c = self.conf
                processed = self.processed

                nt = 2 if low_performance else 0

                with tmpdir.as_cwd():
                        tpr = grompp(f, c, processed, prefix="reference")
                        reference_trr = mdrun(tpr, prefix="reference", nt=nt)
                        df1 = rerun_energy(tpr, reference_trr, prefix="reference", nt=nt)

                        scaled = "scaled.top"
                        kwargs = dict(banned_lines='', topfile=processed, outfile=scaled,
                                      scale_lipids=1.0, scale_protein=1.0)
                        scaling.partial_tempering(**kwargs)

                        assert os.path.exists(scaled), "failed to produce {0}".format(scaled)

                        tpr = grompp(f, c, scaled, prefix="scaled")
                        df2 = rerun_energy(tpr, reference_trr, prefix="scaled", nt=nt)

                        assert_frame_equal(df1, df2, check_names=True, check_like=True)

                        scaled = "scaled.top"
                        kwargs = dict(banned_lines='', topfile=processed,
                                      outfile=scaled, scale_lipids=1.0, scale_protein=0.5)
                        scaling.partial_tempering(**kwargs)
                        tpr = grompp(f, c, scaled, prefix="scaled", maxwarn=1)
                        df3 = rerun_energy(tpr, reference_trr, prefix="scaled", nt=nt)
                        # print(df1, df1.columns)
                        # print(df3, df3.columns)
                        unscaled_terms = ['Time (ps)', 'Improper Dih.']
                        scaled_terms = ['Proper Dih.']

                        assert_frame_equal(df1[unscaled_terms],
                                           df3[unscaled_terms],
                                           check_names=True, check_like=True)
                        assert_frame_equal(df1[scaled_terms],
                                           2*df3[scaled_terms],
                                           check_names=True, check_like=True)



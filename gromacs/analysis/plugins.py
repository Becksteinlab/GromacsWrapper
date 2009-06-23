# $Id$
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""\
Analysis plugins
================

Mixin classes for core.Simulation that provide code to analyze
trajectory data.

See docs in gromacs.analysis.core for preliminary API.

**ALPHA STATUS**
"""
__docformat__ = "restructuredtext en"

import sys
import os.path
import warnings
import subprocess

import gromacs

from core import AttributeDict, Worker, Plugin
import mindist


# Worker classes that are registered via Plugins (see below)
# ----------------------------------------------------------
# These must be defined before the plugins.

class _CysAccessibility(Worker):
    """Analysis of Cysteine accessibility."""

    plugin_name = "CysAccessibility"

    def __init__(self,**kwargs):
        """Set up  customized Cysteine accessibility analysis.

        :Arguments:
        cysteines       list of *all* resids (eg from the sequence) that are used as
                        labels or in the form 'Cys<resid>'. MUST BE PROVIDED.
        cys_cutoff      cutoff in nm for the minimum S-OW distance [1.0]

        Note that *all* Cys residues in the protein are analyzed. Therefore,
        the list of cysteine labels should contain as many entries as there are
        cysteines in the protein.
        """
        super(_CysAccessibility,self).__init__(**kwargs)
        
        # specific setup
        cysteines = kwargs.pop('cysteines',None)     # sequence resids as labels (NOT necessarily Gromacs itp)
        cys_cutoff = kwargs.pop('cys_cutoff', 1.0)   # nm

        # super class do this before doing anything else (maybe not important anymore)
        super(_CysAccessibility,self).__init__(**kwargs)

        self.location = 'accessibility'     # directory under topdir()
        self.results = AttributeDict()
        self.parameters = AttributeDict()

        try:
            self.parameters.cysteines = map(int, cysteines)  # sequence resids
        except (TypeError,ValueError):
            raise ValueError("Keyword argument cysteines MUST be set to sequence of resids.")
        self.parameters.cysteines.sort()                 # sorted because make_ndx returns sorted order
        self.parameters.cutoff = cys_cutoff
        self.parameters.ndx = self.plugindir('cys.ndx')
        # output filenames for g_dist, indexed by Cys resid
        self.parameters.filenames = dict(\
            [(resid, self.plugindir('Cys%d_OW_dist.txt.bz2' % resid))
             for resid in self.parameters.cysteines])
        # default filename for the combined plot
        self.parameters.figname = self.plugindir('mindist_S_OW')

    # override 'API' methods of base class
        
    def run(self,**kwargs):
        return self.run_g_dist_cys(**kwargs)

    def analyze(self,**kwargs):
        return self.analyze_cys()

    # specific methods

    def make_index_cys(self):
        """Make index file for all cysteines and water oxygens. 

        **NO SANITY CHECKS**: The SH atoms are simply labelled consecutively
        with the resids from the cysteines parameter.
        """
        commands_1 = ['keep 0', 'del 0', 'r CYSH & t S', 'splitres 0', 'del 0']  # CYS-S sorted by resid
        commands_2 = ['t OW', 'q']                                               # water oxygens
        commands = commands_1[:]
        for groupid, resid in enumerate(self.parameters.cysteines):
            commands.append('name %(groupid)d Cys%(resid)d'  % vars())           # name CYS-S groups canonically
        commands.extend(commands_2)
        return gromacs.make_ndx(f=self.simulation.tpr, o=self.parameters.ndx, 
                                input=commands, stdout=None)

    def run_g_dist_cys(self,cutoff=None,**gmxargs):
        """Run ``g_dist -dist cutoff`` for each cysteine and save output for further analysis."""

        if cutoff is None:
            cutoff = self.parameters.cutoff
        else:
            self.parameters.cutoff = cutoff    # record cutoff used

        ndx = self.parameters.ndx
        if not os.path.isfile(ndx):
            warnings.warn("Cysteine index file %r missing: running 'make_index_cys'." % ndx)
            self.make_index_cys()

        for resid in self.parameters.cysteines:
            groupname = 'Cys%(resid)d' % vars()
            commands = [groupname, 'OW']
            filename = self.parameters.filenames[resid]
            if self.check_file_exists(filename, resolve='warning'):
                continue
            print "run_g_dist: %(groupname)s --> %(filename)r" % vars()
            sys.stdout.flush()
            datafile = open(filename, 'w')
            try:
                p = gromacs.g_dist.Popen(
                    s=self.simulation.tpr, f=self.simulation.xtc, n=ndx, dist=cutoff, input=commands, 
                    stderr=None, stdout=subprocess.PIPE, **gmxargs)
                compressor = subprocess.Popen(['bzip2', '-c'], stdin=p.stdout, stdout=datafile)
                p.communicate()
            finally:
                datafile.close()

    def analyze_cys(self):
        """Mindist analysis for all cysteines. Returns results for interactive analysis."""        
        results = AttributeDict()
        for resid in self.parameters.cysteines:
            groupname = 'Cys%(resid)d' % vars()    # identifier should be a valid python variable name
            results[groupname] = self._mindist(resid)
        self.results = results
        return results

    def _mindist(self,resid):
        """Analyze minimum distance for resid."""
        filename = self.parameters.filenames[resid]
        return mindist.Mindist(filename,cutoff=self.parameters.cutoff)



# Public classes that register the worker classes
#------------------------------------------------

class CysAccessibility(Plugin):
    """\
    'CysAccessibility' plugin
    =========================
    
    For each frame of a trajectory, the shortest distance of all water oxygens
    to all cysteine sulphur atoms is computed. For computational efficiency,
    only distances smaller than a cutoff are taken into account. A histogram of
    the distances shows how close water molecules can get to cysteines. The
    closest approach distance should be indicative of the reactivity of the SH
    group with crosslinking agents.
    """
    plugin_name = "CysAccessibility"   # XXX: these get overwritten when mixing-in
    plugin_class = _CysAccessibility   # (find a better way to do this..only tested with single mixin yet)



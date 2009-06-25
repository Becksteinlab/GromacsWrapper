# $Id$
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
Distance plugin
===============

Time series of a selected set of distances.

Plugin class
------------

.. autoclass:: Distances
   :members: plugin_name, plugin_class
   :undoc-members:

Worker class
------------

The worker class performs the analysis.

.. autoclass:: _Distances
   :members:


"""
__docformat__ = "restructuredtext en"

import sys
import os.path
import warnings
import subprocess

import gromacs

from gromacs.analysis.core import AttributeDict, Worker, Plugin
import mindist


# Worker classes that are registered via Plugins (see below)
# ----------------------------------------------------------
# These must be defined before the plugins.

class _Distances(Worker):
    """Analysis of distances."""

    plugin_name = "Distances"

    def __init__(self,**kwargs):
        """Set up  customized distance analysis.

        :Arguments:
           A : list
             First group of atoms.
           B : list
             Second group of atoms.
        """
        super(_Distances,self).__init__(**kwargs)
        
        # specific setup
        A = kwargs.pop('A',None)
        B = kwargs.pop('B', None)

        # super class do this before doing anything else (maybe not important anymore)
        super(_Distances,self).__init__(**kwargs)

        self.location = 'distances'     # directory under topdir()
        self.results = AttributeDict()
        self.parameters = AttributeDict()

        self.parameters.A = A
        self.parameters.B = B
        self.parameters.ndx = self.plugindir('distances.ndx')
        # output filenames for g_dist, indexed by Cys resid
        self.parameters.filenames = self.plugindir('dist.txt.bz2')   # ??? needed???

        # default filename for the combined plot
        self.parameters.figname = self.plugindir('distances')

    # override 'API' methods of base class
        
    def run(self,**kwargs):
        return self.run_g_dist(**kwargs)

    def analyze(self,**kwargs):
        return self.analyze_dist()

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

class Distances(Plugin):
    """*Distances* plugin.

    Distances between all members of groups A and B are calculated for
    each time step and written to file
    (:meth:`_Distances.run`). Additional analysis is deferred to the
    :meth:`_Distances.analyze` call.
    """
    plugin_name = "Distances"   # XXX: these get overwritten when mixing-in
    plugin_class = _Distances   # (find a better way to do this..only tested with single mixin yet)



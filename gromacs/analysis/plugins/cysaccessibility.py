# GromacsWrapper plugin: cysaccessibility.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
CysAccessibility plugin
=======================

Cysteine accessibility is analyzed by histogramming the distance of
shortest approach of water molecules to the sulfhydryl group of Cys.

See class docs for more details.

.. note::

   This plugin is the canonical example for how to structure plugins that
   conform to the plugin API (see docs :mod:`gromacs.analysis.core` for
   details).

Plugin class
------------

.. autoclass:: CysAccessibility
   :members: worker_class
   :undoc-members:

Worker class
------------

The worker class performs the analysis.

.. autoclass:: _CysAccessibility
   :members:


"""
from __future__ import with_statement

__docformat__ = "restructuredtext en"

import sys
import os.path
import warnings
import subprocess

import gromacs
from gromacs.utilities import AttributeDict
from gromacs.analysis.core import Worker, Plugin
import dist


# Worker classes that are registered via Plugins (see below)
# ----------------------------------------------------------
# These must be defined before the plugins.

class _CysAccessibility(Worker):
    """Analysis of Cysteine accessibility."""

    def __init__(self,**kwargs):
        """Set up  customized Cysteine accessibility analysis.

        :Arguments:
           cysteines : list
              list of *all* resids (eg from the sequence) that are used as
              labels or in the form 'Cys<resid>'. (**required**)
           cys_cutoff : number
              cutoff in nm for the minimum S-OW distance [1.0]

        Note that *all* Cys residues in the protein are analyzed. Therefore,
        the list of cysteine labels *must* contain as many entries as there are
        cysteines in the protein. There are no sanity checks.
        """
        # specific arguments: take them before calling the super class that
        # does not know what to do with them
        cysteines = kwargs.pop('cysteines',None)     # sequence resids as labels (NOT necessarily Gromacs itp)
        cys_cutoff = kwargs.pop('cys_cutoff', 1.0)   # nm

        # super class init: do this before doing anything else
        # (also sets up self.parameters and self.results)
        super(_CysAccessibility,self).__init__(**kwargs)

        # process specific parameters now
        try:
            self.parameters.cysteines = map(int, cysteines)  # sequence resids
        except (TypeError,ValueError):
            raise ValueError("Keyword argument cysteines MUST be set to sequence of resids.")
        self.parameters.cysteines.sort()                 # sorted because make_ndx returns sorted order
        self.parameters.cutoff = cys_cutoff

        if not self.simulation is None:
            self._register_hook()

    def _register_hook(self, **kwargs):
        """Run when registering; requires simulation."""

        super(_CysAccessibility, self)._register_hook(**kwargs)
        assert not self.simulation is None

        # filename of the index file that we generate for the cysteines
        self.parameters.ndx = self.plugindir('cys.ndx')
        # output filenames for g_dist, indexed by Cys resid
        self.parameters.filenames = {resid: self.plugindir('Cys%d_OW_dist.txt.bz2' % resid)
             for resid in self.parameters.cysteines}
        # default filename for the combined plot
        self.parameters.figname = self.figdir('mindist_S_OW')


    # override 'API' methods of base class
        
    def run(self, cutoff=None, force=False, **gmxargs):
        """Run ``g_dist -dist cutoff`` for each cysteine and save output for further analysis.

        By default existing data files are *not* overwritten. If this
        behaviour is desired then set the *force* = ``True`` keyword
        argument.

        All other parameters are passed on to :func:`gromacs.g_dist`.
        """

        if cutoff is None:
            cutoff = self.parameters.cutoff
        else:
            self.parameters.cutoff = cutoff    # record cutoff used

        ndx = self.parameters.ndx     # could move this into _register_hook?
        if not os.path.isfile(ndx):
            warnings.warn("Cysteine index file %r missing: running 'make_index_cys'." % ndx)
            self.make_index_cys()

        for resid in self.parameters.cysteines:
            groupname = 'Cys%(resid)d' % vars()
            commands = [groupname, 'OW']
            filename = self.parameters.filenames[resid]
            if not force and self.check_file_exists(filename, resolve='warning'):
                continue
            print "run_g_dist: %(groupname)s --> %(filename)r" % vars()
            sys.stdout.flush()
            with open(filename, 'w') as datafile:
                p = gromacs.g_dist.Popen(
                    s=self.simulation.tpr, f=self.simulation.xtc, n=ndx, dist=cutoff, input=commands, 
                    stderr=None, stdout=subprocess.PIPE, **gmxargs)
                compressor = subprocess.Popen(['bzip2', '-c'], stdin=p.stdout, stdout=datafile)
                p.communicate()

    def analyze(self,**kwargs):
        """Mindist analysis for all cysteines. Returns results for interactive analysis."""        

        results = AttributeDict()
        for resid in self.parameters.cysteines:
            groupname = 'Cys%(resid)d' % vars()    # identifier should be a valid python variable name
            results[groupname] = self._mindist(resid)
        self.results = results
        return results

    def plot(self, **kwargs):
        """Plot all results in one graph, labelled by the result keys.

        :Keywords:
           figure
               - ``True``: save figures in the given formats
               - "name.ext": save figure under this filename (``ext`` -> format)
               - ``False``: only show on screen
           formats : sequence
               sequence of all formats that should be saved [('png', 'pdf')]
           plotargs    
               keyword arguments for pylab.plot()
        """

        import pylab
        figure = kwargs.pop('figure', False)
        extensions = kwargs.pop('formats', ('pdf','png'))
        for name,result in self.results.items():
            kwargs['label'] = name
            result.plot(**kwargs)
        pylab.legend(loc='best')
        if figure is True:
            for ext in extensions:
                self.savefig(ext=ext)
        elif figure:
            self.savefig(filename=figure)

    
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

    def _mindist(self,resid):
        """Analyze minimum distance for resid."""
        filename = self.parameters.filenames[resid]
        return dist.Mindist(filename,cutoff=self.parameters.cutoff)


# Public classes that register the worker classes
#------------------------------------------------

class CysAccessibility(Plugin):
    """*CysAccessibility* plugin.
    
    For each frame of a trajectory, the shortest distance of all water oxygens
    to all cysteine sulphur atoms is computed. For computational efficiency,
    only distances smaller than a cutoff are taken into account. A histogram of
    the distances shows how close water molecules can get to cysteines. The
    closest approach distance should be indicative of the reactivity of the SH
    group with crosslinking agents.

    .. class:: CysAccessibility(cysteines, [cys_cutoff[, name[, simulation]]])
    
    :Arguments:
        name : string
            plugin name (used to access it)
        simulation : instance
            The :class:`gromacs.analysis.Simulation` instance that owns the plugin.
        cysteines : list
            list of *all* resids (eg from the sequence) that are used as
            labels or in the form 'Cys<resid>'. (**required**)
        cys_cutoff : number
            cutoff in nm for the minimum S-OW distance [1.0]

    Note that *all* Cys residues in the protein are analyzed. Therefore, the
    list of cysteine labels *must* contain as many entries as there are
    cysteines in the protein. There are no sanity checks.

    """
    worker_class = _CysAccessibility



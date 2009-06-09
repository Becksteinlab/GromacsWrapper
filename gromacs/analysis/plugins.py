# $Id$
"""\
Analysis plugins
================

Mixin classes for core.Simulation that provide code to analyze
trajectory data.

See docs in gromacs.analysis.core for preliminary API.
ALPHA.
"""

import sys
import os.path
import warnings
import subprocess

from core import AttributeDict
import gromacs


# worker classes (used by the plugins)
# (These must be defined before the plugins.)

class _CysAccessibility(object):
    """Analysis of Cysteine accessibility."""

    plugin_name = "CysAccessibility"

    def __init__(self,**kwargs):
        """Set up  customized Cysteine accessibility analysis.

        :Arguments:
        cysteines       list of resids (eg from the sequence) that are used as
                        labels or in the form 'Cys<resid>'. []
        cys_cutoff      cutoff in nm for the minimum S-OW distance [1.0]                        
        """
        # general
        self.simulation = kwargs.pop('simulation',None)  # required (but kw for super & friends)
        assert self.simulation != None
        
        # specific setup
        cysteines = kwargs.pop('cysteines',[])       # sequence resids as labels (NOT necessarily Gromacs itp)
        cys_cutoff = kwargs.pop('cys_cutoff', 1.0)   # nm

        # super class do this before doing anything else (maybe not important anymore)
        super(_CysAccessibility,self).__init__(**kwargs)

        self.location = 'accessibility'     # directory under topdir()
        self.results = AttributeDict()
        self.parameters = AttributeDict()

        self.parameters.cysteines = map(int, cysteines)  # sequence resids
        self.parameters.cysteines.sort()                 # sorted because make_ndx returns sorted order
        self.parameters.cutoff = cys_cutoff
        self.parameters.ndx = self.plugindir('cys.ndx')
        # output filenames for g_dist, indexed by Cys resid
        self.parameters.filenames = dict(\
            [(resid, self.plugindir('Cys%d_OW_dist.txt.bz2' % resid))
             for resid in self.parameters.cysteines])
        # default filename for the combined plot
        self.parameters.figname = self.plugindir('mindist_S_OW')
        
    def topdir(self, *args):
        return self.simulation.topdir(*args)
    
    def plugindir(self, *args):
        return self.topdir(self.location, *args)

    def run(self,**kwargs):
        return self.run_g_dist_cys(**kwargs)

    def analyze(self,**kwargs):
        return self.analyze_cys()

    # specific methods

    def make_index_cys(self):
        """Make index file for all cysteines and water oxygens. 
        NO SANITY CHECKS
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
        return Mindist(filename,cutoff=self.parameters.cutoff)


# plugins:
# registers a worker class in Simulation.plugins and adds a pointer to Simulation to worker

class Plugin(object):
    """Plugin mixin classes are derived from Plugin. 

    They only register the actual plugin in the plugins dictionary.
    """    
    # XXX: gets overwritten with multiple plugin mixins --- do something else!
    plugin_name = None     # name of the plugin
    plugin_class = None    # actual plugin class (typically name with leading underscore)

    def __init__(self,**kwargs):
        """Registers the plugin with the simulation class.
        :Arguments:
        <plugin_name>      a dictionary named like the plugin is taken to include
                           keyword arguments to initialize the plugin
        **kwargs           all other kwargs are passed along                           
        """
        plugin_args = kwargs.pop(self.plugin_name,{})  # must be a dict named like the plugin
        plugin_args['simulation'] = self               # allows access of plugin to globals
        super(Plugin, self).__init__(**kwargs)
        self.plugins[self.plugin_name] = self.plugin_class(**plugin_args)



class CysAccessibility(Plugin):
    plugin_name = "CysAccessibility"   # XXX: these get overwritten when mixing-in
    plugin_class = _CysAccessibility   # (find a better way to do this..only tested with single mixin yet)



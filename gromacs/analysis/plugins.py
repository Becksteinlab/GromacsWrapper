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

class Plugin(object):
    """Plugin mixin classes are derived from Plugin. 

    They only register the actual plugin in the plugins dictionary.
    """    
    plugin_name = None     # name of the plugin
    plugin_class = None    # actual plugin class (typically name with leading underscore)

    def __init__(self,**kwargs):
        """Registers the plugin with the simulation class.
        :Arguments:
        <plugin_name>      a dictionary named like the plugin is taken to include
                           keyword arguments to initialize the plugin
        **kwargs           all other kwargs are passed along                           
        """
        plugin_args = kwargs.pop(plugin_name,None)  # must be a dict named like the plugin
        plugin_args['simulation'] = self            # allows access of plugin to globals
        super(Plugin, self).__init__(**kwargs)
        self.plugins[plugin_name] = plugin_class(**plugin_args)

class XXXCysAccessibility(Plugin):
    plugin_name = "CysAccessibility"
###XXX    plugin_class = _CysAccessibility

# * This is not a well-thought out plan. If the method resolution order is messed
#   up then some data structures (location, parameters, results, _analysis) are
#   not available yet that are needed. 
# * It would be better to encapsulate this in a separate class as to avoid
#   clashes between analysis methods of different plugins.
#
# We'll probably rewrite all this eventually....

#class _CysAccessibility(Plugin):
class CysAccessibility(object):
    """Analysis of Cysteine accessibility. Use as mixin class for Simulation"""

    plugin_name = "CysAccessibility"   # <--- crap, if we mixin more than one Plugin this will be overwritten!!

    def __init__(self,**kwargs):
        """Set up  customized Cysteine accessibility analysis.

        :Arguments:
        cysteines       list of resids (eg from the sequence) that are used as
                        labels or in the form 'Cys<resid>'. []
        cys_cutoff      cutoff in nm for the minimum S-OW distance [1.0]                        
        """
        ## self.simulation = kwargs.pop('simulation',None)  # required

        cysteines = kwargs.pop('cysteines',[])       # sequence resids as labels (NOT necessarily Gromacs itp)
        cys_cutoff = kwargs.pop('cys_cutoff', 1.0)   # nm

        # super class Simulation must come before doing anything else
        super(CysAccessibility,self).__init__(**kwargs)

        # plugin_name = 'CysAccessibility'
        # this should be set by deriving from a base class, currently hard coded for Mhp1
        self.location[self.plugin_name] = 'accessibility'     # directory under topdir()
        self.results[self.plugin_name] = AttributeDict()
        self.parameters[self.plugin_name] = AttributeDict()

        self.parameters[self.plugin_name].cysteines = cysteines
        self.parameters[self.plugin_name].cutoff = cys_cutoff
        self.parameters[self.plugin_name].ndx = self.topdir(self.location[self.plugin_name],'cys.ndx')
        # output filenames for g_dist, indexed by Cys resid
        self.parameters[self.plugin_name].filenames = dict(\
            [(resid, self.plugindir(self.plugin_name, 'Cys%d_OW_dist.txt.bz2' % resid))
             for resid in self.parameters[self.plugin_name].cysteines])
        # default filename for the combined plot
        self.parameters[self.plugin_name].figname = self.plugindir(self.plugin_name, 'mindist_S_OW')

        # analysis method
        self._analysis[self.plugin_name] = self.analyze_cys
        
    def make_index_cys(self):
        """Make index file for all cysteines and water oxygens. 
        NO SANITY CHECKS
        """
        commands_1 = ['keep 0', 'del 0', 'r CYSH & t S', 'splitres 0', 'del 0']
        commands_2 = ['t OW', 'q']
        commands = commands_1[:]
        for groupid, resid in enumerate(self.parameters[self.plugin_name].cysteines):
            commands.append('name %(groupid)d Cys%(resid)d'  % vars())
        commands.extend(commands_2)
        # print "DEBUG: "+'\n'.join(commands)
        return gromacs.make_ndx(f=self.tpr, o=self.parameters[self.plugin_name].ndx, 
                                input=commands, stdout=None)

    def run_g_dist_cys(self,cutoff=None,**gmxargs):
        """Run ``g_dist -dist cutoff`` for each cysteine and save output for further analysis."""

        if cutoff is None:
            cutoff = self.parameters[self.plugin_name].cutoff
        else:
            self.parameters[self.plugin_name].cutoff = cutoff    # record cutoff used

        ndx = self.parameters[self.plugin_name].ndx
        if not os.path.isfile(ndx):
            warnings.warn("Cysteine index file %r missing: running 'make_index_cys'." % ndx)
            self.make_index_cys()

        for resid in self.parameters[self.plugin_name].cysteines:
            groupname = 'Cys%(resid)d' % vars()
            commands = [groupname, 'OW']
            filename = self.parameters[self.plugin_name].filenames[resid]
            print "run_g_dist: %(groupname)s --> %(filename)r" % vars()
            sys.stdout.flush()
            datafile = open(filename, 'w')
            try:
                p = gromacs.g_dist.Popen(
                    s=self.tpr, f=self.xtc, n=ndx, dist=cutoff, input=commands, 
                    stderr=None, stdout=subprocess.PIPE, **gmxargs)
                compressor = subprocess.Popen(['bzip2', '-c'], stdin=p.stdout, stdout=datafile)
                p.communicate()
            finally:
                datafile.close()

    def analyze_cys(self):
        """Mindist analysis for all cysteines. Returns results for interactive analysis."""        
        results = AttributeDict()
        for resid in self.parameters[self.plugin_name].cysteines:
            groupname = 'Cys%(resid)d' % vars()    # identifier should be a valid python variable name
            results[groupname] = self._mindist(resid)
        self.results[self.plugin_name] = results
        return results


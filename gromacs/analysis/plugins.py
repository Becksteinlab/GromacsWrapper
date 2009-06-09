# $Id$
"""\
Analysis plugins
================

Mixin classes for core.Simulation that provide code to analyze
trajectory data.

See docs in gromacs.analysis.core for preliminary API.
ALPHA.
"""

from core import AttributeDict
import gromacs

class Plugin(object):
    """Empty class to nicely label all plugins.

    If we want to add some general code for all plugins in the future it can go
    here.
    """    

# * This is not a well-thought out plan. If the method resolution order is messed
#   up then some data structures (location, parameters, results, _analysis) are
#   not available yet that are needed. 
# * It would be better to encapsulate this in a separate class as to avoid
#   clashes between analysis methods of different plugins.
#
# We'll probably rewrite all this eventually....

class CysAccessibility(Plugin):
    """Analysis of Cysteine accessibility. Use as mixin class for Simulation"""

    plugin_name = "CysAccessibility"

    def __init__(self,**kwargs):
        """Set up  customized Cysteine accessibility analysis.

        :Arguments:
        cysteines       list of resids (eg from the sequence) that are used as
                        labels or in the form 'Cys<resid>'. []
        cys_cutoff      cutoff in nm for the minimum S-OW distance [1.0]                        
        """

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
            [(resid, self.topdir(self.location[self.plugin_name],
                                 'Cys%d_OW_dist.txt.bz2' % resid)) 
             for resid in self.parameters[self.plugin_name].cysteines])

        # analysis method
        self._analysis[self.plugin_name] = self.analyze_cys
        
    def make_cys_index(self):
        """Make index file for all cysteines and water oxygens."""
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
            warnings.warn("Cysteine index file %r missing: running 'make_cys_ndx'." % ndx)
            self.make_cys_index()

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


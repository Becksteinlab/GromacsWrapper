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
import tempfile

import gromacs

from gromacs.analysis.core import AttributeDict, Worker, Plugin
import mindist


# Worker classes that are registered via Plugins (see below)
# ----------------------------------------------------------
# These must be defined before the plugins.

class _Distances(Worker):
    """Analysis of distances.

    First generate index files with the groups of interest::

      from gromacs.cbook import IndexBuilder
      A_ndx = IndexBuilder(tpr, ['@a 62549 & r NA'], names=['Na1'], offset=-9, out_ndx='Na1.ndx', name_all="NA1").combine()
      B_ndx = IndexBuilder('md_posres.pdb', ['S312:OG','T313:OG1','A38:O','I41:O','A309:O'], offset=-9, out_ndx='Na1_site.ndx', name_all="Na1_site").combine()

    
      
      

    """

    plugin_name = "Distances"

    def __init__(self,**kwargs):
        """Set up  customized distance analysis.

        :Arguments:
           A : index group name
             First group of atoms.
           B : index group name
             Second group of atoms.
           ndx : index filename or list
             All index files that contain the A and B groups.
        """
        # specific setup
        A = kwargs.pop('A',None)
        B = kwargs.pop('B', None)
        ndx = kwargs.pop('ndx', None)

        # super class do this before doing anything else (maybe not important anymore)
        super(_Distances,self).__init__(**kwargs)

        self.location = 'distances'     # directory under topdir()
        self.parameters.A = A
        self.parameters.B = B
        self.parameters.ndx = ndx



        # output filenames for g_dist, indexed by Cys resid
        self.parameters.filenames = self.plugindir('dist.txt.bz2')   # ??? needed???

        # default filename for the combined plot
        self.parameters.figname = self.plugindir('distances')

    # override 'API' methods of base class
        
    def run(self,**kwargs):
        return self.run_g_mindist(**kwargs)

    def analyze(self,**kwargs):
        return self.analyze_dist()

    # specific methods

    def run_g_mindist(self,cutoff=None,**gmxargs):
        """Run ``g_mindist ``  for all distances between A and B atoms and save output for further analysis."""

        

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



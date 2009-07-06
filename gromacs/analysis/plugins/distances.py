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
   :members: plugin_name, worker_class
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
from gromacs.utilities import AttributeDict, XVG
from gromacs.analysis.core import Worker, Plugin


# Worker classes that are registered via Plugins (see below)
# ----------------------------------------------------------
# These must be defined before the plugins.

class _Distances(Worker):
    """Analysis of distances.

    See :class:`Distances` for usage.
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
           cutoff : float
             A contact is recorded if the distance is <cutoff [0.6 nm]
        """
        # specific setup
        A = kwargs.pop('A',None)
        B = kwargs.pop('B', None)
        ndx = kwargs.pop('ndx', None)
        cutoff = kwargs.pop('cutoff', 0.6)   # default: 0.6 nm

        # super class: do this before doing anything else
        super(_Distances,self).__init__(**kwargs)

        self.location = 'distances'     # directory under topdir()
        self.parameters.A = A
        self.parameters.B = B
        self.parameters.ndx = ndx
        self.parameters.cutoff = cutoff


        # output filenames for g_dist
        self.parameters.filenames = {'contacts': self.plugindir('contacts.xvg'),
                                     'distance': self.plugindir('distance.xvg'),
                                     }

        # default filename for the combined plot
        self.parameters.figname = self.plugindir('distances')

    # override 'API' methods of base class
        
    def run(self,**kwargs):
        return self.run_g_mindist(**kwargs)

    def analyze(self,**kwargs):
        return self.analyze_dist()

    # specific methods

    def run_g_mindist(self,**gmxargs):
        """Run ``g_mindist `` to compute distances between A and B groups.

        Additional arguments can be provided (e.g. ``-b`` or ``-e``)
        but an error will result if one tries to set parameters that
        are already being set by the method itself such as ``-s`` or
        ``-d``; one must to provide the appropriate values to the
        class constructor.

        If the primary output file already exists then no data are generated
        and the method returns immediately.        
        """
        if self.check_file_exists(self.parameters.filenames['distance'],
                                  resolve='warn'):
            return
            
        indexgroups = [x for x in (self.parameters.A, self.parameters.B)
                       if not x is None]
        gromacs.g_mindist(s=self.simulation.tpr, n=self.parameters.ndx,
                          f=self.simulation.xtc, d=self.parameters.cutoff,
                          od=self.parameters.filenames['distance'],
                          on=self.parameters.filenames['contacts'],
                          o=True, _or=True,    # just generate with default names
                          input=indexgroups,
                          **gmxargs)

    def analyze_dist(self):
        """Make data available as numpy arrays."""        
        results = AttributeDict()

        for name, f in self.parameters.filenames.items():
            results[name] = XVG(f).asarray()

        self.results = results
        return results




# Public classes that register the worker classes
#------------------------------------------------

class Distances(Plugin):
    """*Distances* plugin.

    Distances between the centers of groups A and B and the number of contacts
    are calculated for each time step and written to files.
    
    First generate index files with the groups of interest::

      from gromacs.cbook import IndexBuilder
      A_grp, A_ndx = IndexBuilder(tpr, ['@a 62549 & r NA'], names=['Na1'], offset=-9, 
                                  out_ndx='Na1.ndx', name_all="NA1").combine()
      B_grp, B_ndx = IndexBuilder(tpr, ['S312:OG','T313:OG1','A38:O','I41:O','A309:O'], offset=-9, 
                                  out_ndx='Na1_site.ndx', name_all="Na1_site").combine()

    Use a :class:`gromacs.analysis.Simulation` class with the :class:`Distances` plugin
    and provide the parameters in the ``Distances`` dict::
 
      from gromacs.analysis import Simulation
      S = Simulation(..., plugins=['Distances', ...], 
                     Distances={'A':A_grp, 'B':B_grp, 'ndx': [A_ndx, B_ndx]})

    """
    plugin_name = "Distances"
    worker_class = _Distances



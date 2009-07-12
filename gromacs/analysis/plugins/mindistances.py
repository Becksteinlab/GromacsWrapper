# $Id$
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
MinDistance plugin
===============

Time series of a selected set of shortest distances.

Plugin class
------------

.. autoclass:: MinDistances
   :members: worker_class
   :undoc-members:

Worker class
------------

The worker class performs the analysis.

.. autoclass:: _MinDistances
   :members:


"""
__docformat__ = "restructuredtext en"

import gromacs
from gromacs.utilities import asiterable
from gromacs.analysis.core import Plugin
from distances import _Distances

# Worker classes that are registered via Plugins (see below)
# ----------------------------------------------------------
# These must be defined before the plugins.

class _MinDistances(_Distances):
    """Analysis of (multiple) minimum distances."""

    #: list of results (not used at the moment, see _register_hook())
    names = ["distance", "contacts"]

    default_plot_columns = Ellipsis   # plot everything by default

    def _register_hook(self, **kwargs):
        """Run when registering; requires simulation.

        Defines output files (note that we overwrite the
        parameters.filenames and figname that super might have set).
        """

        super(_MinDistances, self)._register_hook(**kwargs)
        assert not self.simulation is None

        # output filenames for g_mindist
        self.parameters.filenames = {'contacts': self.plugindir('contacts.xvg'),
                                     'distance': self.plugindir('distance.xvg'),
                                     }
        # default filename for the combined plot
        self.parameters.figname = self.figdir('mindist')
    
    def run(self,**kwargs):
        """Run ``g_mindist `` to compute distances between multiple groups.

        Additional arguments can be provided (e.g. ``-b`` or ``-e``)
        but an error will result if one tries to set parameters that
        are already being set by the method itself such as ``-s`` or
        ``-d``; one must to provide the appropriate values to the
        class constructor.

        If the primary output file already exists then no data are generated
        and the method returns immediately unless one sets *force* = ``True``.
        """
        force = kwargs.pop('force',False)
        if not force and \
           self.check_file_exists(self.parameters.filenames['distance'],resolve='warn'):
            return
        indexgroups = self.parameters.indexgroups
        ngroups = len(indexgroups) - 1    # number of secondary groups
        kwargs.setdefault('o', None)     # set to True if default output is required, or
        kwargs.setdefault('_or', None)   #   add filenames to self.parameters.filenames
        gromacs.g_mindist(s=self.simulation.tpr, n=self.parameters.ndx,
                          f=self.simulation.xtc, d=self.parameters.cutoff,
                          od=self.parameters.filenames['distance'],
                          on=self.parameters.filenames['contacts'],
                          ng=ngroups, input=indexgroups,
                          **kwargs)





# Public classes that register the worker classes
#------------------------------------------------

class MinDistances(Plugin):
    """*MinDistances* plugin.

    The minimum distances between the members of at least two index
    groups and the number of contacts are calculated for each time
    step and written to files.

    .. class:: Distances(groups, ndx, [cutoff, [, name[, simulation]]])
    
    :Arguments:
        name : string
            plugin name (used to access it)
        simulation : instance
            The :class:`gromacs.analysis.Simulation` instance that owns the plugin.
        groups : list of index group names
            The first entry is the *primary group*. All other entries
            are *secondary groups* and the plugin calculates the minimum distance
            between members of the primary group and the members of each
            secondary group.
        ndx : index filename or list
            All index files that contain the listed groups.
        cutoff : float
            A contact is recorded if the distance is <cutoff [0.6 nm]

    Example:
    
    Generate index files with the groups of interest, for instance
    with :class:`gromacs.cbook.IndexBuilder`::

      from gromacs.cbook import IndexBuilder
      A_grp, A_ndx = IndexBuilder(tpr, ['@a 62549 & r NA'], names=['Na1_ion'], offset=-9, 
                                  out_ndx='Na1.ndx', name_all="Na1").combine()
      B = IndexBuilder(tpr, ['S312:OG','T313:OG1','A38:O','I41:O','A309:O'], offset=-9, 
                            out_ndx='Na1_site.ndx', name_all="Na1_site")
      B_grp, B_ndx = B.combine()                            
      all_ndx_files = [A_ndx, B_ndx]

    To calculate the distance between "Na1" and the "Na1_site", create an instance with
    the appropriate parameters and add them to a :class:`gromacs.analysis.Simulation` instance::

      dist_Na1_site = Distances(name='Dsite', groups=['Na1', 'Na1_site'], ndx=all_ndx_files)
      S.add_plugin(dist_Na1_site)

    To calculate the individual distances::

      dist_Na1_res = Distances(name='Dres', groups=['Na1']+B.names, ndx=all_ndx_files)
      S.add_plugin(dist_Na1_res)

    (Keeping the second IndexBuilder instance ``B`` allows us to directly
    use all groups without typing them, ``B.names = ['A309_O', 'S312_OG', 'I41_O',
    'T313_OG1', 'A38_O']``.)
    

    """
    worker_class = _MinDistances



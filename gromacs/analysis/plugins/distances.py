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
   :members: worker_class
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
import numpy

import gromacs
from gromacs.utilities import AttributeDict, asiterable
from gromacs.formats import XVG
from gromacs.analysis.core import Worker, Plugin


# Worker classes that are registered via Plugins (see below)
# ----------------------------------------------------------
# These must be defined before the plugins.

class _Distances(Worker):
    """Analysis of distances.

    See :class:`Distances` for usage.

    Also used as a base class for :class:`mindistances._MinDistances`.
    """

    #: list of results (not used at the moment, see _register_hook())
    names = ["distance"]

    #: dict of labels for the plot x-axis; one for each result
    xlabels = {"distance": r"time $t/$ns",
               "contacts": r"time $t/$ns",
               }
    #: dict of labels for the plot y-axis; one for each result
    ylabels = {"distance": r"distance $d/$nm",
               "contacts": r"contacts $N$",
               }

    default_plot_columns = [0, 1]   # plot time and distance only

    def __init__(self,**kwargs):
        """Set up  customized distance analysis.

        :Arguments:
           groups : list of index group names
             The first entry is the *primary group*. All other entries
             are *secondary groups* and the plugin calculates the minimum distance
             between members of the primary group and the members of each
             secondary group.
           ndx : index filename or list
             All index files that contain the listed groups.
           cutoff : float
             A contact is recorded if the distance is <cutoff [0.6 nm]
        """
        # Note: this init is also used for mindistances._MinDistance
        #       thus we have the add. file contacts that is not used otherwise
        # specific setup
        indexgroups = kwargs.pop('groups',None)
        if indexgroups is None or len(indexgroups) < 2 or type(indexgroups) is str:
            raise ValueError("groups must be a list with at least a primary and secondary group")
        ndx = kwargs.pop('ndx', None)
        cutoff = kwargs.pop('cutoff', 0.6)   # default: 0.6 nm

        # super class: do this before setting any instance attributes
        #  sets self.simulation if available!
        super(_Distances,self).__init__(**kwargs)

        self.parameters.indexgroups = indexgroups
        self.parameters.ndx = ndx
        self.parameters.cutoff = cutoff

        if not self.simulation is None:
            self._register_hook()

    def _register_hook(self, **kwargs):
        """Run when registering; requires simulation.

        Defines output files (note that we overwrite the
        parameters.filenames and figname that super might have set).
        """

        super(_Distances, self)._register_hook(**kwargs)
        assert not self.simulation is None

        # output filenames for g_dist
        self.parameters.filenames = {
            'distance': self.plugindir('distance.xvg'),
            }

        # default filename for the combined plot
        self.parameters.figname = self.figdir('distances')


    # override 'API' methods of base class
    
    def run(self,**kwargs):
        """Run ``g_dist `` to compute distances between A and B groups.

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
           self.check_file_exists(self.parameters.filenames['distance'], resolve='warn'):
            return
        indexgroups = self.parameters.indexgroups
        ngroups = len(indexgroups) - 1    # number of secondary groups
        if ngroups != 1:
            raise ValueError("g_dist can only compute the distance between a primary and a secondary group")
        gromacs.g_dist(s=self.simulation.tpr, n=self.parameters.ndx, f=self.simulation.xtc,
                       o=self.parameters.filenames['distance'],
                       input=indexgroups,
                       **kwargs)

    def analyze(self,**kwargs):
        """Make data files available as numpy arrays."""        
        results = AttributeDict()
        for name, f in self.parameters.filenames.items():
            results[name] = XVG(f)
        self.results = results
        return results

    def plot(self, names=None, **kwargs):
        """Plot the selected data.

        :Arguments:
           names : string or list
              Selects which results should be plotted. ``None`` plots all
              in separate graphs.
           columns : list
              Which columns to plot; typically the default is ok.
           figure
               - ``True``: save figures in the given formats
               - "name.ext": save figure under this filename (``ext`` -> format)
               - ``False``: only show on screen
           formats : sequence
               sequence of all formats that should be saved [('png', 'pdf')]
           callbacks : dict
               **hack**: provide a dictionary that contains callback functions
               to customize the plot. They will be called at the end of
               generating a subplot and must be indexed by *name*. They will
               be called with the keyword arguments *name* and *axis*
               (current subplot axis object)::

                    callback(name=name, axis=ax)
           kwargs
              All other keyword arguments are directly passed to 
              meth:`gromacs.formats.XVG.plot`.
        """
        import pylab

        figure = kwargs.pop('figure', False)
        extensions = kwargs.pop('formats', ('pdf','png'))
        callbacks = kwargs.pop('callbacks', None)
        def ps2ns(a):
            """Transform first column (in ps) to ns."""
            _a = numpy.array(a, copy=True)
            _a[0] *= 0.001
            return _a
        kwargs.setdefault('transform', ps2ns)
        kwargs.setdefault('columns', self.default_plot_columns)

        if names is None:
            names = self.results.keys()
        names = asiterable(names)  # this is now a list (hopefully of strings)
        ngraphs = len(names)
        for plotNum, name in enumerate(names):
            plotNum += 1
            ax = pylab.subplot(1, ngraphs, plotNum)
            try:
                data = self.results[name].plot(**kwargs)   # results are XVG objects with plot method
            except KeyError:
                ax.close()
                raise KeyError('name = %r not known, choose one of %r' % (name, self.results.keys()))
            #pylab.title(r'Distances: %s' % name)
            pylab.xlabel(self.xlabels[name])
            pylab.ylabel(self.ylabels[name])
            
            # hack: callbacks for customization
            if not callbacks is None:
                try:
                    callbacks[name](name=name, axis=ax)
                except KeyError:
                    pass

        # pylab.legend(loc='best')
        if figure is True:
            for ext in extensions:
                self.savefig(ext=ext)
        elif figure:
            self.savefig(filename=figure)


                           

# Public classes that register the worker classes
#------------------------------------------------

class Distances(Plugin):
    """*Distances* plugin.

    The distance between the center of mass of two index groups are
    calculated for each time step and written to files.

    .. class:: Distances(groups, ndx, [cutoff, [, name[, simulation]]])
    
    :Arguments:
        name : string
            plugin name (used to access it)
        simulation : instance
            The :class:`gromacs.analysis.Simulation` instance that owns the plugin.
        groups : list of index group names
            The first entry is the *primary group*, the second is the
            *secondary group.
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
    

    """
    worker_class = _Distances



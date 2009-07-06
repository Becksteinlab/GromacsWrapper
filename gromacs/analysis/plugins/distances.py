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
import numpy

import gromacs
from gromacs.utilities import AttributeDict, XVG, asiterable
from gromacs.analysis.core import Worker, Plugin


# Worker classes that are registered via Plugins (see below)
# ----------------------------------------------------------
# These must be defined before the plugins.

class _Distances(Worker):
    """Analysis of distances.

    See :class:`Distances` for usage.
    """

    plugin_name = "Distances"
    #: list of results (not used at the moment, see __init__)
    names = ["distance", "contacts"]
    #: dict of labels for the plot x-axis; one for each result
    xlabels = {"distance": r"time $t/$ns",
               "contacts": r"time $t/$ns",
               }
    #: dict of labels for the plot y-axis; one for each result
    ylabels = {"distance": r"distance $d/$nm",
               "contacts": r"contacts $N$",
               }

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
        kwargs.setdefault('o', True)      # just generate with default names
        kwargs.setdefault('_or', True)    # just generate with default names
        gromacs.g_mindist(s=self.simulation.tpr, n=self.parameters.ndx,
                          f=self.simulation.xtc, d=self.parameters.cutoff,
                          od=self.parameters.filenames['distance'],
                          on=self.parameters.filenames['contacts'],
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
              meth:`gromacs.utilities.XVG.plot`.
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

        if names is None:
            names = self.results.keys()
        names = asiterable(names)  # this is now a list (hopefully of strings)
        ngraphs = len(names)
        for plotNum, name in enumerate(names):
            plotNum += 1
            ax = pylab.subplot(1, ngraphs, plotNum)
            try:
                data = self.results[name].plot(**kwargs)
            except KeyError:
                ax.close()
                raise KeyError('name = %r not known, choose one of %r' % (name, self.results.keys()))
            #pylab.title(r'Distances: %s' % name)
            pylab.xlabel(self.xlabels[name])
            pylab.ylabel(self.ylabels[name])
            
            # hack: callbacks for customization
            if not callbacks is None:
                callbacks[name](name=name, axis=ax)

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



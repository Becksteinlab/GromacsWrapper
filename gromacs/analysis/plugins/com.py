# plugin: com.py
# Copyright (c) 2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
Centre of mass
==============

Calculate the centre of mass of index groups.

Plugin class
------------

.. autoclass:: COM
   :members: worker_class
   :undoc-members:


Worker class
------------

The worker class performs the analysis.

.. autoclass:: _COM
   :members:


"""
from __future__ import with_statement

__docformat__ = "restructuredtext en"

import os.path
import warnings

import numpy

import gromacs
from gromacs.utilities import AttributeDict, asiterable
from gromacs.analysis.core import Worker, Plugin

import logging
logger = logging.getLogger('gromacs.analysis.plugins.com')

# Worker classes that are registered via Plugins (see below)
# ----------------------------------------------------------
# These must be defined before the plugins.

class _COM(Worker):
    """COM worker class."""

    def __init__(self,**kwargs):
        """Set up COM analysis.

        :Keywords:
           *group_names*
               list of index group names
           *ndx*
               index file if groups are not in the default index
           *offset*
               add the *offset* to the residue numbers in *helixtable* [0]
           *name*
               plugin name [COM]
           *simulation*
               The :class:`gromacs.analysis.Simulation` instance that
               owns the plugin [None]
        """
        group_names = asiterable(kwargs.pop('group_names', []))
        ndx = kwargs.pop('ndx', None)
        offset = kwargs.pop('offset', 0)

        super(_COM, self).__init__(**kwargs)
        
        self.parameters.group_names = group_names
        self.parameters.offset = offset
        self.ndx = ndx

        if not self.simulation is None:
            self._register_hook()

    def _register_hook(self, **kwargs):
        """Run when registering; requires simulation."""

        super(_COM, self)._register_hook(**kwargs)
        assert not self.simulation is None

        if self.ndx is None:
            self.ndx = self.simulation.ndx

        self.parameters.filenames = {                     # result xvg files
            'com': self.plugindir('com.xvg'),
            }

        # default filename for the plots -- not used
        self.parameters.fignames = {
            'com': self.figdir('com'),
            }
            
    def run(self, force=None, **gmxargs):
        """Analyze trajectory and write COM file.

        All three components of the COM coordinate are written.

        :Arguments:
          - *force*: ``True`` does analysis and overwrites existing files
          - *gmxargs*: additional keyword arguments for :func:`gromacs.g_bundle` 
        """
        gmxargs['com'] = True
        gmxargs['mol'] = False
        gmxargs['ng'] = len(self.parameters.group_names)
        gmxargs['x'] = True
        gmxargs['y'] = True
        gmxargs['z'] = True

        if gmxargs['ng'] == 0:
            errmsg = "No index group name(s) provided. Use group_name with the constructor."
            logger.error(errmsg)
            raise ValueError(errmsg)

        if self.check_file_exists(self.parameters.filenames['com'], resolve='warning', force=force):
            return

        logger.info("Analyzing COM ...")
        f = self.parameters.filenames
        gromacs.g_traj(s=self.simulation.tpr, f=self.simulation.xtc, n=self.ndx,
                       ox=f['com'],  input=self.parameters.group_names,  **gmxargs)


    def analyze(self,**kwargs):
        """Collect output xvg files as :class:`gromacs.formats.XVG` objects.

        :Returns:  a dictionary of the results and also sets ``self.results``.
        """        
        from gromacs.formats import XVG

        logger.info("Preparing COM graphs as XVG objects.")        
        results = AttributeDict( (k, XVG(fn)) for k,fn in self.parameters.filenames.items() )
        self.results = results
        return results

    def plot(self, **kwargs):
        """Plot all results in one graph, labelled by the result keys.

        :Keywords:
           observables
              select one or more of the stored results. Can be a list
              or a string (a key into the results dict). ``None``
              plots everything [``None``]           
           figure
               - ``True``: save figures in the given formats
               - "name.ext": save figure under this filename (``ext`` -> format)
               - ``False``: only show on screen [``False``]
           formats : sequence
               sequence of all formats that should be saved [('png', 'pdf')]
           plotargs    
               keyword arguments for pylab.plot()
        """

        import pylab
        figure = kwargs.pop('figure', False)
        observables = asiterable(kwargs.pop('observables', self.results.keys()))
        extensions = kwargs.pop('formats', ('pdf','png'))

        for name in observables:
            result = self.results[name]
            try:
                result.plot(**kwargs)      # This requires result classes with a plot() method!!
            except AttributeError:
                warnings.warn("Sorry, plotting of result %(name)r is not implemented" % vars(),
                              category=UserWarning)

        # quick labels -- relies on the proper ordering
        labels = [str(n)+" "+dim for n in self.parameters.group_names
                  for dim in 'xyz']
        if not kwargs.get('columns', None) is None:
            # select labels according to columns; only makes sense
            # if plotting against the time (col 0)
            if kwargs['columns'][0] == 0:
                labels = numpy.array([None]+labels)[kwargs['columns'][1:]]
            else:
                labels = ()

        pylab.legend(labels, loc='best')
        if figure is True:
            for ext in extensions:
                self.savefig(ext=ext)
        elif figure:
            self.savefig(filename=figure)

    


# Public classes that register the worker classes
#------------------------------------------------

class COM(Plugin):
    """*COM* plugin.

    :func:`gromacs.g_bundle` helix analysis

    .. class:: COM([name[, simulation]]])
    
    :Arguments:
       *name*
           plugin name [COM]
       *simulation*
           The :class:`gromacs.analysis.Simulation` instance that owns
           the plugin. [None]
    """
    worker_class = _COM



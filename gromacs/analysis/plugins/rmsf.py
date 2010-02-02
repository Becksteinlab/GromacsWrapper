# $Id$
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
RMSF calculation
================

Calculation of the root mean square fluctuations from a trajectory.

Plugin class
------------

.. autoclass:: RMSF
   :members: worker_class
   :undoc-members:

Worker class
------------

The worker class performs the analysis.

.. autoclass:: _RMSF
   :members:


"""
from __future__ import with_statement

__docformat__ = "restructuredtext en"

import os.path
import warnings

import gromacs
from gromacs.utilities import AttributeDict
from gromacs.analysis.core import Worker, Plugin

import logging
logger = logging.getLogger('gromacs.analysis.plugins.rmsf')

# Worker classes that are registered via Plugins (see below)
# ----------------------------------------------------------
# These must be defined before the plugins.

class _RMSF(Worker):
    """RMSF worker class."""

    def __init__(self,**kwargs):
        """Set up RMSF analysis.

        This is the worker class; this is where all the real analysis is done.

        """
        super(_RMSF, self).__init__(**kwargs)
        if not self.simulation is None:
            self._register_hook()

    def _register_hook(self, **kwargs):
        """Run when registering; requires simulation."""

        super(_RMSF, self)._register_hook(**kwargs)
        assert not self.simulation is None

        self.parameters.filenames = {
            'RMSF': self.plugindir('rmsf.xvg'),
            'RMSD': self.plugindir('rmsdev.xvg'),
            }
        # default filename for the plot
        self.parameters.figname = self.figdir('rmsf')

    def run(self, force=False, group='C-alpha', **gmxargs):
        """Analyze trajectory and write RMSF files.

        The run method typically processes trajectories and writes data files.

        :Arguments:
          - *group*: index group (eg C-alpha or Protein)
          - *force*: do analysis and overwrite existing files
          - *gmxargs*: additional keyword arguments for :func:`gromacs.g_rmsf` (e.g. res=True)
        """
        if not self.check_file_exists(self.parameters.filenames['RMSF'], resolve='warning') or force:
            logger.info("Analyzing RMSF...")
            gromacs.g_rmsf(s=self.simulation.tpr, f=self.simulation.xtc, fit=True, 
                           o=self.parameters.filenames['RMSF'],
                           od=self.parameters.filenames['RMSD'], 
                           input=[group], **gmxargs)

    def analyze(self,**kwargs):
        """Collect output xvg files as :class:`gromacs.formats.XVG` objects.

        :Returns:  a dictionary of the results and also sets ``self.results``.
        """        
        from gromacs.formats import XVG

        logger.info("Preparing RMSF graphs as XVG objects.")
        results = AttributeDict(RMSF=XVG(self.parameters.filenames['RMSF']),
                                RMSD=XVG(self.parameters.filenames['RMSD']))
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
            try:
                result.plot(**kwargs)      # This requires result classes with a plot() method!!
            except AttributeError:
                warnings.warn("Sorry, plotting of result %(name)r is not implemented" % vars(),
                              category=UserWarning)                
        pylab.legend(loc='best')
        if figure is True:
            for ext in extensions:
                self.savefig(ext=ext)
        elif figure:
            self.savefig(filename=figure)

    


# Public classes that register the worker classes
#------------------------------------------------

class RMSF(Plugin):
    """*RMSF* plugin.
    
    Compute the root mean square fluctuations (RMSF) of the C-alpha
    atoms. The trajectory is always fitted to the reference structure
    in the tpr file.

    .. class:: RMSF([name[, simulation]])
    
    :Arguments:
        *name* : string
            plugin name (used to access it)
        *simulation* : instance
            The :class:`gromacs.analysis.Simulation` instance that owns the plugin.

    """
    worker_class = _RMSF



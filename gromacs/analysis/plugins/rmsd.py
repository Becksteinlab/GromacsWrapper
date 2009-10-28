# $Id$
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
RMSD calculation
================

Calculation of the root mean square distance of a protein structure
over the course of a MD simulation.


Plugin class
------------

.. autoclass:: RMSD
   :members: worker_class
   :undoc-members:

Worker class
------------

The worker class performs the analysis.

.. autoclass:: _RMSD
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
logger = logging.getLogger('gromacs.analysis.plugins.rmsd')

# Worker classes that are registered via Plugins (see below)
# ----------------------------------------------------------
# These must be defined before the plugins.

class _RMSD(Worker):
    """RMSD worker class."""

    def __init__(self,**kwargs):
        """Set up RMSD analysis.

        This is the worker class; this is where all the real analysis is done.
        """
        super(_RMSD, self).__init__(**kwargs)
        if not self.simulation is None:
            self._register_hook()

    def _register_hook(self, **kwargs):
        """Run when registering; requires simulation."""

        super(_RMSD, self)._register_hook(**kwargs)
        assert not self.simulation is None

        self.parameters.filenames = {
            'RMSD': self.plugindir('rmsd.xvg'),
            }
        # default filename for the plot
        self.parameters.figname = self.figdir('rmsd')

    def run(self, force=False, group='C-alpha', **gmxargs):
        """Analyze trajectory and write RMSD files.

        Each frame is fit to the structure in the tpr, based on the
        C-alpha atoms. The RMSd is alculated for the supplied index
        group *group*.

        :Arguments:
          - *group*: index group for RMSD calculation (eg C-alpha or Protein)
          - *force*: do analysis and overwrite existing files
          - *gmxargs*: additional keyword arguments for :func:`gromacs.g_rms` (e.g. res=True)
        """
        if not self.check_file_exists(self.parameters.filenames['RMSD'], resolve='warning') or force:
            logger.info("Analyzing RMSD...")
            gromacs.g_rms(s=self.simulation.tpr, f=self.simulation.xtc, fit="rot+trans", 
                          o=self.parameters.filenames['RMSD'],
                          input=['C-alpha', group], **gmxargs)

    def analyze(self,**kwargs):
        """Collect output xvg files as :class:`gromacs.formats.XVG` objects.

        :Returns:  a dictionary of the results and also sets ``self.results``.
        """        
        from gromacs.formats import XVG

        logger.info("Preparing RMSD graphs as XVG objects.")
        results = AttributeDict(RMSD=XVG(self.parameters.filenames['RMSD']))
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

class RMSD(Plugin):
    """*RMSD* plugin.
    
    Calculation of the root mean square distance (RMSD) of a protein
    structure over the course of a MD simulation.

    The trajectory is always fitted to the reference structure in the
    tpr file.

    .. class:: RMSD([name[, simulation]])
    
    :Arguments:
        *name* : string
            plugin name (used to access it)
        *simulation* : instance
            The :class:`gromacs.analysis.Simulation` instance that owns the plugin.

    """
    worker_class = _RMSD



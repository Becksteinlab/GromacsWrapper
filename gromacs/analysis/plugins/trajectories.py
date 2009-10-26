# $Id$
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
Trajectories
============

Write centered trajectories. Trajectories will be stored under the
base dir / trj.  

- full fitxy 
- fitxy at 100 ps intervals


Plugin class
------------

.. autoclass:: Trajectories
   :members: worker_class
   :undoc-members:

Worker class
------------

The worker class performs the analysis.

.. autoclass:: _Trajectories
   :members:


"""
from __future__ import with_statement

__docformat__ = "restructuredtext en"

import os
import errno
import warnings

import gromacs
from gromacs.utilities import AttributeDict
from gromacs.analysis.core import Worker, Plugin


# Worker classes that are registered via Plugins (see below)
# ----------------------------------------------------------
# These must be defined before the plugins.

class _Trajectories(Worker):
    """Trajctories worker class."""

    def __init__(self,**kwargs):
        """Set up  Trajectories

        :Arguments:
           None at the moment, everything is hard coded.
        """
        # specific arguments: take them before calling the super class that
        # does not know what to do with them
        dt = kwargs.pop('dt', 100)   # reduced xtc: write steps every dt ps

        # super class init: do this before doing anything else
        # (also sets up self.parameters and self.results)
        super(_Trajectories, self).__init__(**kwargs)

        # process specific parameters now and set instance variables
        self.parameters.dt = dt

        # self.simulation might have been set by the super class
        # already; just leave this snippet at the end. Do all
        # initialization that requires the simulation class in the
        # _register_hook() method.
        if not self.simulation is None:
            self._register_hook()

    def _register_hook(self, **kwargs):
        """Run when registering; requires simulation."""

        super(_Trajectories, self)._register_hook(**kwargs)
        assert not self.simulation is None        

        xtcdir,xtcname = os.path.split(self.simulation.xtc)
        xtcbasename, xtcext = os.path.splitext(xtcname)
        trjdir = os.path.join(xtcdir, 'trj')
        fitxy_gro = os.path.join(trjdir, xtcbasename+'_fitxy'+'.gro')
        fitxy_pdb = os.path.join(trjdir, xtcbasename+'_fitxy'+'.pdb')
        fitxy_xtc = os.path.join(trjdir, xtcbasename+'_fitxy'+xtcext)
        fitxydt_xtc = os.path.join(trjdir, (xtcbasename+'_fitxy_dt%dps'+xtcext) % self.parameters.dt)

        # make sure trjdir exists
        try:
            os.makedirs(trjdir)
        except OSError, err:
            if err.errno == errno.EEXIST:
                pass

        self.parameters.trjdir = trjdir        
        self.parameters.filenames = {
            'gro': fitxy_gro,
            'pdb': fitxy_pdb,
            'fitxy': fitxy_xtc,
            'fitxydt': fitxydt_xtc,
            }
        self.parameters.ndx = self.plugindir('nowater.ndx')
        self.parameters.groups = {'nowater': 'nowater'}

    # override 'API' methods of base class
        
    def run(self, force=False, **gmxargs):
        """Write new trajectories"""

        logger.info("Writing fitted GRO file...")
        gromacs.cbook.trj_fitandcenter(xy=True, s=self.simulation.tpr, f=self.simulations.xtc, 
                                       o=self.parameters.filenames['gro'], dump=0)

        logger.info("Writing fitted PDB file...")
        gromacs.cbook.trj_fitandcenter(xy=True, s=self.simulation.tpr, f=self.simulations.xtc, 
                                       o=self.parameters.filenames['pdb'], dump=0)

        logger.info("Writing fitted xtc file (frame every %d ps)..." % self.parameters.dt)
        gromacs.cbook.trj_fitandcenter(xy=True, s=self.simulation.tpr, f=self.simulations.xtc, 
                                       o=self.parameters.filenames['fitxydt'],
                                       dt=self.parameters.dt)

        logger.info("Writing fitted xtc file (all frames)...")
        gromacs.cbook.trj_fitandcenter(xy=True, s=self.simulation.tpr, f=self.simulations.xtc, 
                                       o=self.parameters.filenames['fitxy'])

        logger.info("New trajectories can be found in %r." % self.parameters.trjdir)

    def analyze(self,**kwargs):
        """No postprocessing."""        
        pass


    def plot(self, **kwargs):
        """No plotting."""
        pass


# Public classes that register the worker classes
#------------------------------------------------

class Trajectories(Plugin):
    """*Trajectories* plugin.

    Write new xy-fitted trajectories (see :func:`gromacs.cbook.trj_fitandcenter`),

    .. class:: Trajectories([name[, simulation]])
    
    :Arguments:
        *name* : string
            plugin name (used to access it)
        *simulation* : instance
            The :class:`gromacs.analysis.Simulation` instance that owns the plugin.

    """
    worker_class = _Trajectories



# plugin for GromacsWrapper: stripwater.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
StripWater
==========

Write a trajectory with all water removed. This uses
:meth:`gromacs.cbook.Transformer.strip_water`.

Plugin class
------------

.. autoclass:: Stripwater
   :members: worker_class
   :undoc-members:

Worker class
------------

The worker class performs the analysis.

.. autoclass:: _StripWater
   :members:


"""
from __future__ import with_statement

__docformat__ = "restructuredtext en"

import os.path
import warnings

import gromacs
import gromacs.cbook
from gromacs.utilities import AttributeDict
from gromacs.analysis.core import Worker, Plugin


# Worker classes that are registered via Plugins (see below)
# ----------------------------------------------------------
# These must be defined before the plugins.

class _StripWater(Worker):
    """StripWater worker class."""

    def __init__(self,**kwargs):
        """Set up  StripWater

        :Arguments:
           *fit*
               one of "xy", "all", ``None``
        """
        # specific arguments: take them before calling the super class that
        # does not know what to do with them
        _fitvalues = ("xy", "all", None)
        parameters = {}
        parameters['fit'] = kwargs.pop('fit',None)            # fitting algorithm
        if not parameters['fit'] in _fitvalues:
            raise ValueError("StripWater: *fit* must be one of %(_fitvalues)r, not %(fit)r." % vars())
        parameters['compact'] = kwargs.pop('compact', False)  # compact+centered ?
        parameters['resn'] = kwargs.pop('resn', 'SOL')        # residue name to be stripped

        # super class init: do this before doing anything else
        # (also sets up self.parameters and self.results)
        super(_StripWater, self).__init__(**kwargs)

        # self.parameters is set up by the base Worker class...
        self.parameters.filenames = AttributeDict()
        self.parameters.update(parameters)

        # self.simulation might have been set by the super class
        # already; just leave this snippet at the end. Do all
        # initialization that requires the simulation class in the
        # _register_hook() method.
        if not self.simulation is None:
            self._register_hook()

    def _register_hook(self, **kwargs):
        """Run when registering; requires simulation."""

        super(_StripWater, self)._register_hook(**kwargs)
        assert not self.simulation is None

        self.transformer = gromacs.cbook.Transformer(
            s=self.simulation.tpr, f=self.simulation.xtc, n=self.simulation.ndx)

    # override 'API' methods of base class
        
    def run(self, force=False, **kwargs):
        """Write new trajectory with water index group stripped.

        kwargs are passed to :meth:`gromacs.cbook.Transformer.strip_water`.

        .. Note:: If set, *dt* is only applied to a fit step; the
                  no-water trajectory is always generated for all time
                  steps of the input.
        """
        dt = kwargs.pop('dt', None)
                
        kwargs.setdefault('compact', self.parameters.compact)
        kwargs.setdefault('resn', self.parameters.resn)

        newfiles = self.transformer.strip_water(**kwargs)
        self.parameters.filenames.update(newfiles)

        if self.parameters.fit != None:
            if self.parameters.fit == "xy": xy = True
            else: xy = False
            
            transformer_nowater = self.transformer.nowater.values()[0]
            transformer_nowater.fit(xy=xy, dt=dt)

    def analyze(self,**kwargs):
        """No postprocessing."""        
        pass


    def plot(self, **kwargs):
        """No plotting."""
        pass


# Public classes that register the worker classes
#------------------------------------------------

class StripWater(Plugin):
    """*StripWater* plugin.

    Write a new trajectory which has the water index group removed. 

    .. class:: StripWater([selection[, name[, simulation]]])
    
    :Arguments:
        *selection*
            optional selection for the water instead of "SOL"
        *name* : string
            plugin name (used to access it)
        *simulation* : instance
            The :class:`gromacs.analysis.Simulation` instance that owns the plugin.

    """
    worker_class = _StripWater



# $Id$
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
StripWater
==========

Write a trajectory with all water removed.


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
from gromacs.utilities import AttributeDict
from gromacs.analysis.core import Worker, Plugin


# Worker classes that are registered via Plugins (see below)
# ----------------------------------------------------------
# These must be defined before the plugins.

class _StripWater(Worker):
    """TEMPLATE worker class."""

    def __init__(self,**kwargs):
        """Set up  StripWater


        :Arguments:
           *fit*
               one of "xy", "all", ``None``
        """
        # specific arguments: take them before calling the super class that
        # does not know what to do with them
        _fitvalues = ("xy", "all", None)
        if not fit in _fitvalues:
            raise ValueError("*fit* must be one of %(_fitvalues)r, not %(fit)r." % vars())
        self.fit = kwargs.pop('fit',None)

        # super class init: do this before doing anything else
        # (also sets up self.parameters and self.results)
        super(_StripWater, self).__init__(**kwargs)

        # process specific parameters now and set instance variables
        # ....

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

        xtcdir,xtcname = os.path.split(self.simulation.xtc)
        xtcbasename, xtcext = os.path.splitext(xtcname)
        newxtc = os.path.join(xtcdir, xtcbasename+'_nowater'+xtcext)
        
        self.parameters.filenames = {'newxtc': newxtc}
        self.parameters.ndx = self.plugindir('nowater.ndx')
        self.parameters.groups = {'nowater': 'nowater'}

    # override 'API' methods of base class
        
    def run(self, force=False, **gmxargs):
        """Write new trajectory with water index group stripped.
        """
        # make ndx without water
        # TODO: ! SOL hard coded

        B = gromacs.cbook.IndexBuilder(struct=self.simulation.tpr, selections=['@! "SOL"'], 
                                       names=[self.groups['nowater']], out_ndx=self.parameters.ndx)
        B.combine()
        if self.fit == None:
            # just strip water
            gromacs.trjconv(s=self.simulation.tpr, f=self.simulation.xtc, n=self.parameters.ndx,
                            o=self.parameters.filenames['newxtc'], input=[self.parameters.groups])
        else:
            # fit and center, then strip water
            xy = {'xy': True, 'all': False}
            gromacs.cbook.trj_fitandcenter(
                xy=xy[self.fit],
                s=self.simulation.tpr, f=self.simulation.xtc, o=self.parameters.filenames['newxtc'],
                n1=self.parameters.ndx, input1=['protein', self.parameters.groups],
                n=None, input=['backbone', 'protein', 'system'])

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



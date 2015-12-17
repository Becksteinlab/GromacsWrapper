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
from gromacs.utilities import AttributeDict, asiterable
from gromacs.analysis.core import Worker, Plugin


# Worker classes that are registered via Plugins (see below)
# ----------------------------------------------------------
# These must be defined before the plugins.

class _StripWater(Worker):
    """StripWater worker class."""

    def __init__(self,**kwargs):
        """Set up  StripWater

        :Arguments:

          *force*
             ``True`` will always regenerate trajectories even if they
             already exist, ``False`` raises an exception, ``None``
             does the sensible thing in most cases (i.e. notify and
             then move on).
          *dt* : float or list of floats
             only write every dt timestep (in ps); if a list of floats is
             supplied, write multiple trajectories, one for each dt.
          *compact* : bool
             write a compact representation
          *centergroup*
             Index group to center on ["Protein"]
          *fit*
             Create an additional trajectory from the stripped one in which
             the *fitgroup* group is rms-fitted to the initial structure. See
             :meth:`gromacs.cbook.Transformer.fit` for details. Useful
             values:

             - "xy" : perform a rot+trans fit in the x-y plane
             - "all": rot+trans
             - ``None``: no fitting

             If *fit* is not supplied then the constructor-default is used
             (:attr:`_StripWater.parameters.fit`).
          *fitgroup*
             Index group to fit to with the *fit* option; must be changed if
             molecule is not a protein and automatically recognized. Also
             consider supplying a custom index file. ["backbone"]
          *resn*
             name of the residues that are stripped (typically it is
             safe to leave this at the default 'SOL')
          *outdir*
             place generated files in *outdir* instead of the same directory
             where the input tpr/xtc lived [``None``]

        .. Note::

           If set, *dt* is only applied to a fit step; the no-water
           trajectory is always generated for all time steps of the
           input.

        """
        # specific arguments: take them before calling the super class that
        # does not know what to do with them
        _fitvalues = ("xy", "all", None)
        parameters = {}
        parameters['fit'] = kwargs.pop('fit',None)            # fitting algorithm
        if not parameters['fit'] in _fitvalues:
            raise ValueError("StripWater: *fit* must be one of %(_fitvalues)r, not %(fit)r." % vars())
	parameters['fitgroup'] = kwargs.pop('fitgroup', "backbone")
	parameters['centergroup'] = kwargs.pop('centergroup', "Protein")
        parameters['compact'] = kwargs.pop('compact', False)  # compact+centered ?
        parameters['resn'] = kwargs.pop('resn', 'SOL')        # residue name to be stripped
        parameters['dt'] = kwargs.pop('dt', None)
        parameters['force'] = kwargs.pop('force', None)
        parameters['outdir'] = kwargs.pop('outdir', None)

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
        if self.simulation is not None:
            self._register_hook()

    def _register_hook(self, **kwargs):
        """Run when registering; requires simulation."""

        super(_StripWater, self)._register_hook(**kwargs)
        assert self.simulation is not None

        trjdir = os.path.dirname(self.simulation.tpr)
        self.transformer = gromacs.cbook.Transformer(dirname=trjdir, outdir=self.parameters.outdir,
            s=self.simulation.tpr, f=self.simulation.xtc, n=self.simulation.ndx)

    # override 'API' methods of base class

    def run(self, **kwargs):
        """Write new trajectory with water index group stripped.

        kwargs are passed to
        :meth:`gromacs.cbook.Transformer.strip_water`. Important
        parameters:

        :Keywords:

          *force*
             ``True`` will always regenerate trajectories even if they
             already exist, ``False`` raises an exception, ``None``
             does the sensible thing in most cases (i.e. notify and
             then move on).
          *dt* : float or list of floats
             only write every dt timestep (in ps); if a list of floats is
             supplied, write multiple trajectories, one for each dt.
          *compact* : bool
             write a compact and centered representation
          *centergroup*
             Index group to center on ["Protein"]
          *fit*
             Create an additional trajectory from the stripped one in which
             the *fitgroup* group is rms-fitted to the initial structure. See
             :meth:`gromacs.cbook.Transformer.fit` for details. Useful
             values:

             - "xy" : perform a rot+trans fit in the x-y plane
             - "all": rot+trans
             - ``None``: no fitting

             If *fit* is not supplied then the constructor-default is used
             (:attr:`_StripWater.parameters.fit`).
          *fitgroup*
             Index group to fit to with the *fit* option; must be changed if
             molecule is not a protein and automatically recognized. Also
             consider supplying a custom index file. ["backbone" or constructor
             supplied]
         *resn*
             name of the residues that are stripped (typically it is
             safe to leave this at the default 'SOL')

        .. Note::

           If set, *dt* is only applied to a fit step; the no-water
           trajectory is always generated for all time steps of the
           input.

        """
        dt = kwargs.pop('dt', self.parameters.dt)
        fit = kwargs.pop('fit', self.parameters.fit)
        fitgroup = kwargs.pop('fitgroup', self.parameters.fitgroup)

        kwargs.setdefault('centergroup', self.parameters.centergroup)
        kwargs.setdefault('compact', self.parameters.compact)
        kwargs.setdefault('resn', self.parameters.resn)
        kwargs.setdefault('force', self.parameters.force)

        newfiles = self.transformer.strip_water(**kwargs)
        self.parameters.filenames.update(newfiles)

        if fit is not None:
            if self.parameters.fit == "xy":
                xy = True
            else:
                xy = False
            transformer_nowater = self.transformer.nowater.values()[0]
            for delta_t in asiterable(dt):
                transformer_nowater.fit(xy=xy, dt=delta_t, fitgroup=fitgroup, force=kwargs['force'])

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

    .. class:: StripWater([resn[,force[,dt[,compact[,fit[,name[,simulation]]]]]]])

    :Arguments:
        *resn*
           name of the residues that are stripped (typically it is
           safe to leave this at the default 'SOL')
        *force*
           ``True`` will always regenerate trajectories even if they
           already exist, ``False`` raises an exception, ``None``
           does the sensible thing in most cases (i.e. notify and
           then move on).
        *dt* : float or list of floats
           only write every dt timestep (in ps); if a list of floats is
           supplied, write multiple trajectories, one for each dt.
        *compact* : bool
           write a compact representation
        *fit*
           Create an additional trajectory from the stripped one in which
           the *fitgroup* group is rms-fitted to the initial structure. See
           :meth:`gromacs.cbook.Transformer.fit` for details. Useful
           values:

             - "xy" : perform a rot+trans fit in the x-y plane
             - "all": rot+trans
             - ``None``: no fitting

           If *fit* is not supplied then the constructor-default is used
           (:attr:`_StripWater.parameters.fit`).
        *fitgroup*
           Index group to fit to with the *fit* option; must be changed if
           molecule is not a protein and automatically recognized. Also
           consider supplying a custom index file. ["backbone"]
        *name* : string
            plugin name (used to access it)
        *simulation* : instance
            The :class:`gromacs.analysis.Simulation` instance that owns the plugin.

    .. Note:: If set, *dt* is only applied to a fit step; the
              no-water trajectory is always generated for all time
              steps of the input.

              The defaults work for typical proteins but for e.g. nucleic
              acids you will need to supply group names in *input* and 
              *fitgroup* and possibly also a custom index file.

    """
    worker_class = _StripWater



# plugin for GromacsWrapper: stripwater.py
# Copyright (c) 2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
ProteinOnly
===========

Write a trajectory with only the protein retained, using
:meth:`gromacs.cbook.Transformer.keep_protein_only`.

Plugin class
------------

.. autoclass:: Stripwater
   :members: worker_class
   :undoc-members:

Worker class
------------

The worker class performs the analysis.

.. autoclass:: _ProteinOnly
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

class _ProteinOnly(Worker):
    """ProteinOnly worker class."""

    def __init__(self,**kwargs):
        """Set up  ProteinOnly

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
          *fit*
             Create an additional trajectory from the stripped one in which
             the Protein group is rms-fitted to the initial structure. See
             :meth:`gromacs.cbook.Transformer.fit` for details. Useful
             values:

             - "xy" : perform a rot+trans fit in the x-y plane
             - "all": rot+trans
             - ``None``: no fitting

             If *fit* is not supplied then the constructore-default is used
             (:attr:`_ProteinOnly.parameters.fit`).
          *keepalso*
             List of literal ``make_ndx`` selections that select additional
             groups of atoms that should also be kept in addition to the
             protein. For example *keepalso* = ['"POPC"', 'resname DRUG'].

        """
        # specific arguments: take them before calling the super class that
        # does not know what to do with them
        _fitvalues = ("xy", "all", None)
        parameters = {}
        parameters['fit'] = kwargs.pop('fit',None)            # fitting algorithm
        if not parameters['fit'] in _fitvalues:
            raise ValueError("ProteinOnly: *fit* must be one of %(_fitvalues)r, not %(fit)r." % vars())
        parameters['compact'] = kwargs.pop('compact', False)  # compact+centered ?
        parameters['dt'] = kwargs.pop('dt', None)
        parameters['force'] = kwargs.pop('force', None)
        parameters['keepalso'] = kwargs.pop('keepalso', None)

        # super class init: do this before doing anything else
        # (also sets up self.parameters and self.results)
        super(_ProteinOnly, self).__init__(**kwargs)

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

        super(_ProteinOnly, self)._register_hook(**kwargs)
        assert not self.simulation is None

        trjdir = os.path.dirname(self.simulation.tpr)
        self.transformer = gromacs.cbook.Transformer(dirname=trjdir,
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
             write a compact representation
          *fit*
             Create an additional trajectory from the stripped one in which
             the Protein group is rms-fitted to the initial structure. See
             :meth:`gromacs.cbook.Transformer.fit` for details. Useful
             values:
               - "xy" : perform a rot+trans fit in the x-y plane
               - "all": rot+trans
               - ``None``: no fitting
             If *fit* is not supplied then the constructore-default is used
             (:attr:`_ProteinOnly.parameters.fit`).
          *keepalso*
              List of ``make_ndx`` selections that should also be kept.

        .. Note::

           If set, *dt* is only applied to a fit step; the no-water
           trajectory is always generated for all time steps of the
           input.
        """
        dt = kwargs.pop('dt', self.parameters.dt)
        fit = kwargs.pop('fit', self.parameters.fit)

        kwargs.setdefault('compact', self.parameters.compact)
        kwargs.setdefault('force', self.parameters.force)
        kwargs.setdefault('keepalso', self.parameters.keepalso)

        newfiles = self.transformer.keep_protein_only(**kwargs)
        self.parameters.filenames.update(newfiles)

        if fit is not None:
            if self.parameters.fit == "xy":
                xy = True
            else:
                xy = False
            transformer_proteinonly = self.transformer.proteinonly.values()[0]
            for delta_t in asiterable(dt):
                transformer_proteinonly.fit(xy=xy, dt=delta_t, force=kwargs['force'])

    def analyze(self,**kwargs):
        """No postprocessing."""
        pass

    def plot(self, **kwargs):
        """No plotting."""
        pass


# Public classes that register the worker classes
#------------------------------------------------

class ProteinOnly(Plugin):
    """*ProteinOnly* plugin.

    Write a new trajectory which has the water index group removed.

    .. class:: ProteinOnly([selection[, name[, simulation[, ...]]]])

    :Arguments:

        *selection*
            optional selection for the water instead of "SOL"
        *name* : string
            plugin name (used to access it)
        *simulation* : instance
            The :class:`gromacs.analysis.Simulation` instance that owns the plugin.
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
           the Protein group is rms-fitted to the initial structure. See
           :meth:`gromacs.cbook.Transformer.fit` for details. Useful
           values:
             - "xy" : perform a rot+trans fit in the x-y plane
             - "all": rot+trans
             - ``None``: no fitting
           If *fit* is not supplied then the constructore-default is used
           (:attr:`_ProteinOnly.parameters.fit`).
        *keepalso*
           List of literal ``make_ndx`` selections that select additional
           groups of atoms that should also be kept in addition to the
           protein. For example *keepalso* = ['"POPC"', 'resname DRUG'].

    """
    worker_class = _ProteinOnly



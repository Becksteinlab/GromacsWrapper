# plugin for GromacsWrapper: stripwater.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
FitCompact
==========

Write a trajectory fitted and centered on the protein, with the box
represented in a "compact" manner.. This uses
:func:`gromacs.cbook.trj_fitandcenter`.

Plugin class
------------

.. autoclass:: FitCompact
   :members: worker_class
   :undoc-members:

Worker class
------------

The worker class performs the analysis.

.. autoclass:: _FitCompact
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

import logging
logger = logging.getLogger('gromacs.analysis.plugins.fitcompact')

# Worker classes that are registered via Plugins (see below)
# ----------------------------------------------------------
# These must be defined before the plugins.

class _FitCompact(Worker):
    """FitCompact worker class."""

    def __init__(self,**kwargs):
        """Set up  FitCompact

        :Arguments:
          *xy* : bool
              Rot+trans fit in 3D (``False``) or trans+rot in x-y plane (``True``) [``False``]
        """
        xy = kwargs.pop('xy', False)
        # maybe also allow input?

        # super class init: do this before doing anything else
        # (also sets up self.parameters and self.results)
        super(_FitCompact, self).__init__(**kwargs)

        # process specific parameters now and set instance variables
        self.parameters.xy = xy

        # self.simulation might have been set by the super class
        # already; just leave this snippet at the end. Do all
        # initialization that requires the simulation class in the
        # _register_hook() method.
        if not self.simulation is None:
            self._register_hook()

    def _register_hook(self, **kwargs):
        """Run when registering; requires simulation."""

        super(_FitCompact, self)._register_hook(**kwargs)
        assert not self.simulation is None

        xtcdir,xtcname = os.path.split(self.simulation.xtc)
        xtcbasename, xtcext = os.path.splitext(xtcname)
        trjdir = os.path.join(xtcdir, 'trj')
        fitcompact_gro = os.path.join(trjdir, xtcbasename+'_fitcompact'+'.gro')
        fitcompact_pdb = os.path.join(trjdir, xtcbasename+'_fitcompact'+'.pdb')
        fitcompact_xtc = os.path.join(trjdir, xtcbasename+'_fitcompact'+xtcext)

        # make sure trjdir exists
        try:
            os.makedirs(trjdir)
        except OSError, err:
            if err.errno == errno.EEXIST:
                pass

        self.parameters.trjdir = trjdir
        self.parameters.filenames = {
            'gro': fitcompact_gro,
            'pdb': fitcompact_pdb,
            'fitcompact': fitcompact_xtc,
            }

    # override 'API' methods of base class

    def run(self, force=False, **gmxargs):
        """Write new trajectories"""

        filename = self.parameters.filenames['gro']
        if not self.check_file_exists(filename, resolve='warning') or force:
            logger.info("Writing fitted GRO file...")
            gromacs.cbook.trj_fitandcenter(xy=self.parameters.xy, s=self.simulation.tpr,
                                           f=self.simulation.xtc,
                                           o=filename, dump=0)

        filename = self.parameters.filenames['pdb']
        if not self.check_file_exists(filename, resolve='warning') or force:
            logger.info("Writing fitted PDB file...")
            gromacs.cbook.trj_fitandcenter(xy=self.parameters.xy, s=self.simulation.tpr,
                                           f=self.simulation.xtc,
                                           o=filename, dump=0)

        filename = self.parameters.filenames['fitcompact']
        if not self.check_file_exists(filename, resolve='warning') or force:
            logger.info("Writing fitted xtc file (all frames)...")
            gromacs.cbook.trj_fitandcenter(xy=self.parameters.xy, s=self.simulation.tpr,
                                           f=self.simulation.xtc,
                                           o=filename)

        logger.info("New trajectories can be found in %r." % self.parameters.trjdir)

    def analyze(self,**kwargs):
        """No postprocessing."""
        pass


    def plot(self, **kwargs):
        """No plotting."""
        pass


# Public classes that register the worker classes
#------------------------------------------------

class FitCompact(Plugin):
    """*FitCompact* plugin.

    Write new trajectories (see :func:`gromacs.cbook.trj_fitandcenter`) where
    the solvent is centered around the protein in a compact representation and
    the system is RMSD-fitted to the protein backbone.

    The output file names are the input names with "_fitcompact" added. For
    instance, if "./MD/md.xtc" is the input trajectory then the output
    trajectory will be "./MD/md_fitcompact.xtc". In addition to the trajectory,
    compact fitted GRO and PDB files will also be produced.

    .. class:: FitCompact([xy[,name[,simulation]]])

    The plugin uses the defaults of
    :func:`gromacs.cbook.trj_fitancenter`; in particular it fits on the
    *backbone*, centers on *protein*, and writes *system*.

    :Arguments:
        *xy* : bool
            Rot+trans fit in 3D (``False``) or trans+rot in x-y plane (``True``) [``False``]
        *name* : string
            plugin name (used to access it)
        *simulation* : instance
            The :class:`gromacs.analysis.Simulation` instance that owns the plugin.

    """
    worker_class = _FitCompact



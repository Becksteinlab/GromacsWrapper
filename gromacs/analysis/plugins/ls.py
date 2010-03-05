# GromacsWrapper plugin: ls.py
# Copyright (c) 2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
ls plugin
=========

This simply lists the files on disk. It is useful for testing the
plugin architecture.


Plugin class
------------

.. autoclass:: Ls
   :members: worker_class
   :undoc-members:

Worker class
------------

The worker class performs the analysis.

.. autoclass:: _Ls
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
logger = logging.getLogger('gromacs.analysis.plugins.ls')

# Worker classes that are registered via Plugins (see below)
# ----------------------------------------------------------
# These must be defined before the plugins.

class _Ls(Worker):
    """ls worker class."""

    def __init__(self,**kwargs):
        """
        :Arguments:
          *kwargs*
             same as ls ??
        """
        # super class init: do this before doing anything else
        # (also sets up self.parameters and self.results)
        super(_Ls, self).__init__(**kwargs)

        # process specific parameters now and set instance variables
        # ....
        # self.parameters.filenames = { 'xxx': 'yyy', ....}
        # ....

        # self.simulation might have been set by the super class
        # already; just leave this snippet at the end. Do all
        # initialization that requires the simulation class in the
        # _register_hook() method.
        if not self.simulation is None:
            self._register_hook()

    def _register_hook(self, **kwargs):
        """Run when registering; requires simulation."""

        super(_Ls, self)._register_hook(**kwargs)
        assert not self.simulation is None


    # override 'API' methods of base class
        
    def run(self, *args, **kwargs):
        """List the contents of the simulation directory.
        """

        from subprocess import call
        lscmd = ['ls', '-la'] + list(args)
        cmd = lscmd + [self.simulation.tpr, self.simulation.xtc]
        with rulify("TPR and XTC"):
            rc = call(cmd)   # just print to screen

        adir = self.simulation.analysis_dir
        cmd = lscmd + [adir]
        with rulify("Analysis dir %(adir)s" % vars()):
            rc = call(cmd)   # just print to screen
        
    def analyze(self,**kwargs):
        pass

    def plot(self, **kwargs):
        pass

from contextlib import contextmanager

@contextmanager
def rulify(header, ncol=79):
    toprule = ncol * '='
    midrule = ncol * '-'
    botrule = toprule
    print toprule
    print header
    print midrule
    try:
        yield None
    finally:
        print botrule
        print

# Public classes that register the worker classes
#------------------------------------------------

class Ls(Plugin):
    """*ls* plugin.
    
    This simply lists the files on disk. It is useful for testing the
    plugin architecture.

    .. class:: Ls([name[, simulation]])
    
    :Arguments:
        *name* : string
            plugin name (used to access it)
        *simulation* : instance
            The :class:`gromacs.analysis.Simulation` instance that owns the plugin.

    """
    worker_class = _Ls



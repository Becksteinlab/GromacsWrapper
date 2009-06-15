# $Id$
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""Core classes for simulation of Gromacs trajectories.

Programming API for plugins
===========================
* Additional analysis capabilities are added with mixin classes (from plugins).
* Derive class for the simulation of interest along the lines of ::

  from gromacs.analysis import Simulation
  from gromacs.analysis.plugins import CysAccessibility

  class Mhp1(Simulation, CysAccessibility):
    def __init__(self,**kwargs):
        kwargs['CysAccessibility'] = {'cysteines': [69, 234, 327]}
        super(Mhp1,self).__init__(**kwargs)

  S = Mhp1(tpr=..., xtc=..., analysisdir=...)
  S.set_default_plugin('CysAccessibility')
  S.run()
  S.analyze()
  S.plot(figure=True)


Plugins
-------

Analysis capabilities can be added by mixing in additional plugins into the
simulation base class. Each plugin registers itself and provides at a minimum
run(), analyze(), and plot() methods.

The plugin class is derived from Plugin and bears the name that is used to
access it. When its __init__ method executes it adds the actual worker class
(typically named with the underscore-prepended name) to the Simulation.plugins
dictionary.

Variables for initializing a plugin a given to the class constructor as a
keyword argument that is named like the plugin and contains a dictionary that
is used as the keyword parameters for the plugin's init.

A plugin **must** obtain a pointer to the Simulation class as the keyword
argument ``simulation`` in order to be able to access simulation-global
parameters such as top directories or input files.

See CysAccessibility and _CysAccessibility as examples.

"""


import sys, os
import errno
import subprocess
import warnings

import gromacs
from mindist import Mindist

class Simulation(object):
    """Class that represents one simulation.

    Analysis capailities are added with mixin classes. 

    NOTE: Only keyword arguments are allowed so that init resolution
    works as expected.
    """
    def __init__(self, **kwargs):
        """Set up a Simulation object; analysis is performed via methods.

        :Arguments:
        tpr             Gromacs tpr file (required)
        xtc             Gromacs trajectory (required)
        ndx             Gromacs index file
        analysisdir     directory under which derived data are stored;
                        defaults to the directory containing the tpr [None]        
        """
        # required files
        self.tpr = kwargs.pop('tpr',None)
        self.xtc = kwargs.pop('xtc',None)
        for v in ('tpr', 'xtc'):            
            self.check_file(v, self.__getattribute__(v))

        self.ndx = kwargs.pop('ndx',None)
        self.analysis_dir = kwargs.pop('analysisdir', os.path.dirname(self.tpr))

        # registry for plugins: This dict is central.
        self.plugins = AttributeDict()   # nicer for interactive use than dict
        self.default_plugin_name = None

        # important: Do not forget to call mixin classes:
        #            each plugin class registers itself in self.plugins[].
        super(Simulation,self).__init__(**kwargs)

        # XXX: does this work (i.e. do we end up here AFTER mixin inits??)
        # convenience: if only a single plugin was registered we default to that one
        if len(self.plugins) == 1:
            self.set_default_plugin(self.plugins.keys()[0])

    def topdir(self,*args):
        """Returns path under self.analysis_dir. Parent dirs are created if necessary."""
        p = os.path.join(self.analysis_dir, *args)
        parent = os.path.dirname(p)
        try:
            os.makedirs(parent)
        except OSError,err:
            if err.errno != errno.EEXIST:
                raise
        return p

    def plugindir(self,plugin_name,*args):
        return self.select_plugin(plugin_name).plugindir(*args)

    def check_file(self,filetype, path):
        if path is None or not os.path.isfile(path):
            raise ValueError("Missing required file %(filetype)r, got %(path)r." % vars())
        return True

    def set_default_plugin(self,plugin_name):
        """Set the plugin that should be used by default.

        If no plugin_name is supplied to run(), analyze() etc. then
        this will be used.
        """
        if plugin_name == None:
            self.default_plugin_name = None
        else:
            self.check_plugin_name(plugin_name)
            self.default_plugin_name = plugin_name
        return self.default_plugin_name

    def has_plugin(self,plugin_name):
        """Returns True if plugin_name is registered."""
        return plugin_name in self.plugins

    def check_plugin_name(self,plugin_name):
        """Raises a ValueError if plugin_name is not registered."""
        if not (plugin_name is None or self.has_plugin(plugin_name)):
            raise ValueError('plugin_name must be None or one of\n%r\n' % self.plugins.keys())

    def select_plugin(self,plugin_name=None):
        """Return valid plugin or the default for plugin_name=None."""
        self.check_plugin_name(plugin_name)
        if plugin_name is None:
            if self.default_plugin_name is None:
                raise ValueError('No default plugin was set.')
            plugin_name = self.default_plugin_name
        return self.plugins[plugin_name]

    def run(self,plugin_name=None,**kwargs):
        """Generate data files as prerequisite to analysis."""
        return self.select_plugin(plugin_name).run(**kwargs)

    def analyze(self,plugin_name=None,**kwargs):
        """Run analysis for the plugin."""
        return self.select_plugin(plugin_name).analyze(**kwargs)    

    def plot(self,plugin_name=None,figure=False,**plotargs):
        """Plot all data for the selected plugin.
        
        plot(plugin_name, **kwargs)

        :Arguments:
        plugin_name   name of the plugin to plot data from
        figure        True: plot to file with default name.
                      string: use this filename (+extension for format)
                      False: only display
        **plotargs    arguments for pylab.plot
        """
        return self.select_plugin(plugin_name).plot(figure=figure,**plotargs)    

    def __str__(self):
        return 'Simulation(tpr=%(tpr)r,xtc=%(xtc)r,analysisdir=%(analysis_dir)r)' % vars(self)
    def __repr__(self):
        return str(self)



class AttributeDict(dict):
    """A dictionary with pythonic access to keys as attributes --- useful for interactive work."""
    def __getattribute__(self,x):
        try:
            return super(AttributeDict,self).__getattribute__(x)
        except AttributeError:
            return self[x]
    def __setattr__(self,name,value):
        try:
            super(AttributeDict,self).__setitem__(name, value)
        except KeyError:
            super(AttributeDict,self).__setattr__(name, value)
            
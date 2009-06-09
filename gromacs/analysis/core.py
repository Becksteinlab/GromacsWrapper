#!/usr/bin/env python
# $Id$

"""Core classes for simulation of Gromacs trajectories.

Programming API for plugins
===========================
* Additional analysis capabilities are added with mixin classes (from plugins).
* Derive class for the simulation of interest along the lines of ::

  from plugins import CysAccessibility
  class Mhp1(Simulation, CysAccessibility):
    def __init__(self,**kwargs):
        kwargs['CysAccessibility'] = {'cysteines': [69, 234, 327]}
        super(Mhp1,self).__init__(**kwargs)

  S = Mhp1(tpr=..., xtc=..., analysisdir=...)
  S.run('CysAccessibility')
  S.analyze('CysAccessibility')
  S.plot('CysAccessibility')


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
import gromacs
from mindist import Mindist
import subprocess
import warnings

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

        # registry for plugins
        self.plugins = AttributeDict()   # nicer for interactive use than dict

        # important: do not forget to call mixin classes
        super(Simulation,self).__init__(**kwargs)


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
        return self.plugins[plugin_name].plugindir(*args)

    def check_file(self,filetype, path):
        if path is None or not os.path.isfile(path):
            raise ValueError("Missing required file %(filetype)r, got %(path)r." % vars())
        return True
    
    def run(self,plugin_name,**kwargs):
        """Generate data files as prerequisite to analysis."""
        return self.plugins[plugin_name].run(**kwargs)

    def analyze(self,plugin_name,**kwargs):
        """Run analysis for the plugin."""
        return self.plugins[plugin_name].analyze(**kwargs)    

    def plot(self,plugin_name,figure=True,**plotargs):
        """Plot all data for the selected plugin.
        
        plot(plugin_name, **kwargs)

        :Arguments:
        plugin_name   name of the plugin to plot data from
        figure        True: plot to file with default name.
                      string: use this filename (+extension for format)
                      False: only display
        **plotargs    arguments for pylab.plot
        """
        # XXX: move plot functionality also into plugins
        import pylab
        for name,result in self.plugins[plugin_name].results.items():
            plotargs['label'] = name
            result.plot(**plotargs)
        pylab.legend(loc='best')
        if figure is True:
            for ext in ('pdf','png'):
                fn = self.plugins[plugin_name].parameters.figname + '.' + ext
                pylab.savefig(fn)
        elif figure:
            pylab.savefig(figure)

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
            

#!/usr/bin/env python
# $Id$

"""Core classes for simulation of Gromacs trajectories.

Programming API for plugins
===========================
* Additional analysis capabilities are added with mixin classes (from plugins).
* Derive class for the simulation of interest along the lines of ::

  class Mhp1(Simulation, CysAccessibility):
    def __init__(self,**kwargs):
        kwargs['cysteines'] = [69, 234, 327]
        super(Mhp1,self).__init__(**kwargs)

  S = Mhp1(tpr=..., xtc=..., analysisdir=...)


Plugins
-------

A plugin has the plugin_name attribute; typically it equald the class name.

Each plugin is responsible for filling in the appropriate dictionaries
that are provided by the Simulation main class. See CysAccessibility as an example.

  # Class dictionaries: each plugin records parameters, results, and location here.
  # Key == plugin_name
  # results: should be objects that have a plot() method
  self.location = AttributeDict()
  self.parameters =  AttributeDict()        
  self.results = AttributeDict()
  self._analysis = AttributeDict()     # record the analysis method; signature call(**kwargs)

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
        # XXX should make different analysis plugins/classes with each in a separate dir
        # XXX sort of done... not sure if the idea with mixin classes is good, though
        # XXX Mabe better to make them individual classes (but how to communicate?)

        # required files
        self.tpr = kwargs.pop('tpr',None)
        self.xtc = kwargs.pop('xtc',None)
        for v in ('tpr', 'xtc'):            
            self.check_file(v, self.__getattribute__(v))

        self.ndx = kwargs.pop('ndx',None)
        self.analysis_dir = kwargs.pop('analysisdir', os.path.dirname(self.tpr))

        # registry for plugins
        self.plugins = {}

        # Class dictionaries: each plugin records parameters, results, and location here.
        # Key == plugin_name
        # results: should be objects that have a plot() method
        self.location = AttributeDict()
        self.parameters =  AttributeDict()        
        self.results = AttributeDict()
        self._analysis = AttributeDict()     # record the analysis method; signature call(**kwargs)

        # important: do not forget to call mixin classes
        super(Simulation,self).__init__(**kwargs)


    def topdir(self,*args):
        """Returns path under self.analysis_dir. Parent dirs are created if necessary."""
        p = os.path.join(self.analysis_dir, *args)
        parent = os.path.dirname(p)
        try:
            os.makedirs(parent)
        except OSError,err:
            if err.errno == errno.EEXIST:
                pass
            else:
                raise
        return p

    def plugindir(self,plugin_name,*args):
        return self.topdir(self.location[self.plugin_name], *args)

    def check_file(self,filetype, path):
        if path is None or not os.path.isfile(path):
            raise ValueError("Missing required file %(filetype)r, got %(path)r." % vars())
        return True
    
    def analyze(self,plugin_name,**kwargs):
        """Run analysis for the plugin."""
        return self._analysis[plugin_name](**kwargs)

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
        import pylab
        for name,result in self.results[plugin_name].items():
            plotargs['label'] = name
            result.plot(**plotargs)
        pylab.legend(loc='best')
        if figure is True:
            for ext in ('pdf','png'):
                fn = self.parameters[plugin_name].figname + '.' + ext
                pylab.savefig(fn)
        elif figure:
            pylab.savefig(figure)


    # semi-generic analysis snippets.... could be cleaned up
    def _mindist(self,resid,plugin_name='CysAccessibility'):
        """Analyze minimum distance for resid."""
        filename = self.parameters[plugin_name].filenames[resid]
        return Mindist(filename,cutoff=self.parameters[plugin_name].cutoff)


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
            

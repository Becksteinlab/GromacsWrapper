# $Id$
"""
``analysis`` -- Analysis Package Overview
=========================================

The :mod:`gromacs.analysis` package is a framework for analyzing Gromacs MD
trajectories. The basic object is the :class:`Simulation` class. For a
particular project one has to derive a class from :class:`Simulation` and add
analysis plugin classes (from :mod:`gromacs.analysis.plugins`) for specific
analysis tasks. This is slightly cumbersome but flexible.

New analysis plugins should follow the API sketched out in
:mod:`gromacs.analysis.core`; see an example for use there.

Right now the number of plugins is limited and simply demonstrates how to use
the framework in principle. If you would like to contribute your own plugins
feel free to send then to the `package author`_. If they have been written
according to the API they will be added to the distribution and of course you
will be acknowledged in the list of plugin authors in
:mod:`gromacs.analysis.plugins`.

.. _`package author`: oliver.beckstein@bioch.ox.ac.uk


Simulation class
----------------

The :class:`Simulation` class is central for doing analysis. The user can
derive a custom analysis class that pre-defines values for plugins as seen in
the `Example`_.

.. autoclass:: Simulation
   :members: __init__, add_plugin, set_default_plugin, run, analyze, plot

Example
-------

Here we analyze *MyProtein*, which has three Cysteines (C96, C243, C372). The
only analysis plugin we want to use is :class:`plugins.CysAccessibility` but
because we might want to do additional analysis in the future we also add the
:class:`plugins.Distances` plugin ::

  **TODO**
  **MUST BE UPDATED**

  from gromacs.analysis import Simulation
  from gromacs.analysis.plugins import CysAccessibility

  class MyProtein(Simulation):
    def __init__(self,**kwargs):
        kwargs['CysAccessibility'] = {'cysteines': [96, 243, 372]}  # pre-sets for CysAccessibility
        super(MyProtein,self).__init__(**kwargs)   # should ALWAYS come last

  S = MyProtein(tpr=..., xtc=..., analysisdir=..., plugins=['CysAccessibility', 'Distances'])
  S.set_plugin('CysAccessibility')          # do CysAccessibility analysis
  S.run()                                   # analyze trajectory and write files
  S.analyze()                               # analyze output files
  S.plot(figure=True)                       # plot and save the figure

Note how the :class:`plugins.CysAccessibility` plugin is initialized via ::

  kwargs['CysAccessibility'] = {'cysteines': [96, 243, 372]}

Other plugins might require no or a very different initialization. See the
plugin documentation for what is required.
"""

__docformat__ = "restructuredtext en"
__all__ = ['Simulation', 'plugins']

from core import Simulation
import plugins

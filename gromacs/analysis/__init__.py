# $Id$
"""
``analysis`` -- Analysis Package Overview
=========================================

The ``analysis`` package is a framework for analyzing Gromacs MD
trajectories. The basic object is the Simulation class. For a particular
project one has to derive a class from simulation and mix-in all analysis
plugin classes that are required. This is slightly cumbersome but flexible. 

New analysis plugins should follow the API sketched out in
``analysis.core``; see an example for use there.

Right now the number of plugins is limited. Feel free to contribute your own by
sending it to the `package author`_. You will be acknowledged in the list of
plugin authors in ``plugin.__init__``.

.. _`package author`: oliver.beckstein@bioch.ox.ac.uk


Simulation class
----------------

The :class:`Simulation` class is central for doing analysis. The user
will have to derive a custom analysis class that mixes
:class:`Simulation` and any plugin classes.

.. autoclass:: Simulation


Example
-------

Here we analyze MyProtein, which has three Cysteines (C96, C243, C372). The
only analysis plugin used is CysAccessibility::

  from gromacs.analysis import Simulation
  from gromacs.analysis.plugins import CysAccessibility

  class MyProtein(Simulation, CysAccessibility):
    def __init__(self,**kwargs):
        kwargs['CysAccessibility'] = {'cysteines': [96, 243, 372]}
        super(MyProtein,self).__init__(**kwargs)

  S = MyProtein(tpr=..., xtc=..., analysisdir=...)
  S.set_default_plugin('CysAccessibility')
  S.run()
  S.analyze()
  S.plot(figure=True)


"""
__docformat__ = "restructuredtext en"
__all__ = ['Simulation', 'plugins']

from core import Simulation
import plugins

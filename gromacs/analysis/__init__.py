# GromacsWrappe analysis module
# Copyright (c) 2009-2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License version 3 (or higher)

"""
:mod:`gromacs.analysis` -- Analysis Package Overview
====================================================

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
   :members: add_plugin, set_plugin, run, analyze, plot

Example
-------

Here we analyze a protein, which has three Cysteines (C96, C243, C372). We
will use the :class:`plugins.CysAccessibility` and the
:class:`plugins.Distances` plugin (arguments for ``Distances`` omitted)::

  from gromacs.analysis import Simulation
  from gromacs.analysis.plugins import CysAccessibility, Distances

  S = Simulation(tpr=..., xtc=..., analysisdir=...,
                 plugins=[('CysAccessibility', {'cysteines': [96, 243, 372]}),
                          Distances(...),
                          ])
  S.set_plugin('CysAccessibility')          # do CysAccessibility analysis
  S.run()                                   # analyze trajectory and write files
  S.analyze()                               # analyze output files
  S.plot(figure=True)                       # plot and save the figure

.. Note:: Absolute paths for the files and *analysisdir* are a good
          idea because many plugins change directory freely.

The plugins can be supplied when the ``Simulation`` object is
constructed, or they can be later added, e.g. ::

  S.add_plugin(Distances(name='Dist2', ...))

This second ``Distances`` analysis would be available with ::

  S.set_plugin('Dist2')

Other plugins might require no or a very different initialization. See the
plugin documentation for what is required.
"""

__docformat__ = "restructuredtext en"
__all__ = ['Simulation', 'plugins']

from core import Simulation
import plugins

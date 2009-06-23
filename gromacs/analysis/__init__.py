# $Id$
"""
Analysis Modules
================

The ``analysis`` package is a framework for analyzing Gromacs MD
trajectories. The basic object is the Simulation class. For a particular
project one has to derive a class from simulation and mix-in all analysis
plugin classes that are required. This is slightly cumbersome but flexible. 

New analysis plugins should follow the API sketched out in
``analysis.core``; see an example for use there.

Right now the number of plugins is limited. Feel free to contribute your own by
sending it to the `package author`_. You will be acknowledged in the list below.

.. _`package author`: oliver.beckstein@bioch.ox.ac.uk


Authors
-------

================       ============================== 
plugin                 author
================       ==============================
CysAccessibility       Oliver Beckstein
================       ==============================

"""
__docformat__ = "restructuredtext en"

__all__ = ['plugins']

from core import Simulation

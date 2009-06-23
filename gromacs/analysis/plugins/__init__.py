# $Id$
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.
"""
Plugin Modules
================

Mixin classes for core.Simulation that provide code to analyze
trajectory data.

See docs in gromacs.analysis.core for preliminary API.

**ALPHA STATUS**


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
__all__ = ['CysAccessibility']

# the plugin can/should mask the package of the same name
from CysAccessibility import CysAccessibility

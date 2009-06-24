# $Id$
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.
"""
``analysis.plugins`` -- Plugin Modules
======================================

Mixin classes for ``core.Simulation`` that provide code to analyze
trajectory data.

New analysis plugins should follow the API sketched out in
``analysis.core``; see an example for use there.

Right now the number of plugins is limited. Feel free to contribute your own by
sending it to the `package author`_. You will be acknowledged in the list below.

.. _`package author`: oliver.beckstein@bioch.ox.ac.uk

Authors
-------

====================  ============================== 
plugin                 author
====================  ==============================
CysAccessibility       Oliver Beckstein
====================  ==============================

Warning   
-------
This is **ALPHA STATUS**. Use at your own risk.
"""
__docformat__ = "restructuredtext en"
__all__ = ['CysAccessibility']

# the plugin can/should mask the package of the same name
from CysAccessibility import CysAccessibility

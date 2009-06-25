# $Id$
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.
"""
:mod:`analysis.plugins` -- Plugin Modules
=========================================

Mixin classes for :class:`gromacs.analysis.core.Simulation` that
provide code to analyze trajectory data.

New analysis plugins should follow the API sketched out in
:mod:`gromacs.analysis.core`; see an example for use there.



List of plugins
---------------

Right now the number of plugins is limited. Feel free to contribute your own by
sending it to the `package author`_. You will be acknowledged in the list below.

.. _`package author`: oliver.beckstein@bioch.ox.ac.uk

====================  =======================  ===============================
plugin                author                   description
====================  =======================  ===============================
CysAccessibility      Oliver Beckstein [#OB]_  estimate accessibility of Cys
                                               residues by water
====================  =======================  ===============================


.. rubric:: Footnotes
.. [#OB] oliver.beckstein@bioch.ox.ac.uk

"""
__docformat__ = "restructuredtext en"
__all__ = ['CysAccessibility']

# the plugin can/should mask the package of the same name
from CysAccessibility import CysAccessibility

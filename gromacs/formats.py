# GromacsWrapper: formats.py
# Copyright (c) 2009-2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
:mod:`gromacs.formats` -- Accessing various files
=================================================

This module contains classes that represent data files on
disk. Typically one creates an instance and

- reads from a file using a :meth:`read` method, or

- populates the instance (in the simplest case with a :meth:`set`
  method) and the uses the :meth:`write` method to write the data to
  disk in the appropriate format.

For function data there typically also exists a :meth:`plot` method
which produces a graph (using matplotlib).

The module defines some classes that are used in other modules; they
do *not* make use of :mod:`gromacs.tools` or :mod:`gromacs.cbook` and
can be safely imported at any time.


Classes
-------

.. autoclass:: XVG
   :members:
.. autoclass:: NDX
   :members:
.. autoclass:: uniqueNDX
   :members:
.. autoclass:: MDP
   :members:
.. autoclass:: ITP
   :members:
.. autoclass:: XPM
   :members:

"""
from __future__ import absolute_import

__docformat__ = "restructuredtext en"

__all__ = ["XVG", "MDP", "NDX", "uniqueNDX", "ITP", "XPM", "TOP"]

from .fileformats import XVG, MDP, NDX, uniqueNDX, ITP, XPM, TOP



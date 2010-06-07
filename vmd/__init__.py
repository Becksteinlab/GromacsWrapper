# $Id: __init__.py 1212 2007-10-24 20:23:03Z oliver$
# vmd remote control --- client/server scripts to run VMD from python
# Copyright (c) 2007-2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Lesser Public License, version 3 or later.
# See COPYING and COPYING.LESSER.
"""
=============
 VMD control
=============

Simple client to transmit Tcl commands to a server running in `VMD`_.

VMD and the server run locally and can be started from the module. Once the
server is running, one can use :class:`vmd.client` to communicate with the server
process via a local socket.

Example
-------

Start a VMD server and connect::

  from vmd.control import *
  VMD = server()
  VMD.command('molecule new load 1AKE')

or start an interactive `Tcl`_ session connected to a running VMD 
server process::

  interactive(host)
  asyncore.loop()      # necessary

See `VMD Tcl Text Commands`_ for all available commands.

.. _VMD: http://www.ks.uiuc.edu/Research/vmd/
.. _Tcl: http://www.tcl.tk/man/
.. _VMD Tcl Text Commands: http://www.ks.uiuc.edu/Research/vmd/current/ug/node107.html

"""
__docformat__ = "restructuredtext en"
__all__ = [ 'control', ]

import control
from control import command as cmd


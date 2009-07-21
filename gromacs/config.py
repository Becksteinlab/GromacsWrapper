# $Id$
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""Configuration for GromacsWrapper.

This is a stub for handling user-settable configuration options. See
pylab's ``rc`` mechanism.
"""

#: preliminary: list the variables in gromacs.tools that should be
#: loaded. Possible values: *gmx_tools*, *gmx_extra_tools*.
load_tools = ['gmx_tools', ]

# XXX: install script via setup.py ... this is a hack:
# part of a python egg
# see http://peak.telecommunity.com/DevCenter/PythonEggs#accessing-package-resources
import os
from pkg_resources import resource_filename
GridMAT_MD = resource_filename(__name__,'external/GridMAT-MD_v1.0.2/GridMAT-MD.pl')
os.chmod(GridMAT_MD, 0755)


#: 3rd party analysis scripts and tools; triplets of 
#:   (script name/path, command name, doc string)
load_scripts = [
    (GridMAT_MD,
     'GridMAT_MD',
     """GridMAT-MD: A Grid-based Membrane Analysis Tool for use with Molecular Dynamics.

*This* ``GridMAT-MD`` is a patched version of ``GridMAT-MD.pl``
v1.0.2, written by WJ Allen, JA Lemkul and DR Bevan. The original
version is available from the `GridMAT-MD`_ home page,
   
.. _`GridMAT-MD`: http://www.bevanlab.biochem.vt.edu/GridMAT-MD/index.html

Please cite 

  W. J. Allen, J. A. Lemkul, and D. R. Bevan. (2009) "GridMAT-MD: A
  Grid-based Membrane Analysis Tool for Use With Molecular Dynamics."
  J. Comput. Chem. 30 (12): 1952-1958.

when using this programme.

Usage:

.. class:: GridMAT_MD(config[,structure])

:Arguments:
- *config* : See the original documentation for a description for the
  configuration file.
- *structure* : A gro or pdb file that overrides the value for
  *bilayer* in the configuration file.
"""),
    ]



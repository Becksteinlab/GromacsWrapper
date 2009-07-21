# $Id$
"""Configuration for GromacsWrapper (not really used yet).

This is a stub for handling user-settable configuration options. See
pylab's ``rc`` mechanism. 
"""

#: preliminary: list the variables in gromacs.tools that should be
#: loaded. Possible values: *gmx_tools*, *gmx_extra_tools*.
load_tools = ['gmx_tools', ]

#: 3rd party analysis scripts and tools; triplets of 
#:   (script name/path, command name, doc string)
load_scripts = [
    ('GridMAT-MDx.pl',
     'GridMAT_MDx',
     """GridMAT-MD: A Grid-based Membrane Analysis Tool for use with Molecular Dynamics.

*GridMAT-MDx* is a patched version of ``GridMAT-MD.pl`` , written by WJ
Allen, JA Lemkul and DR Bevan. It is available from the `GridMAT-MD`_
home page,
   
.. _`GridMAT-MD`: http://www.bevanlab.biochem.vt.edu/GridMAT-MD/index.html

The patch can be obtained from Oliver Beckstein <orbeckst@gmail.com>;
it allows overriding the structure file name from the command line.

Please cite 

  W. J. Allen, J. A. Lemkul, and D. R. Bevan. (2009) "GridMAT-MD: A
  Grid-based Membrane Analysis Tool for Use With Molecular Dynamics."
  J. Comput. Chem. 30 (12): 1952-1958.

when using this programme.

Usage:

.. class:: GridMAT_MDx(config[,structure])

:Arguments:
- *config* : See the original documentation for a description for the
  configuration file.
- *structure* : A gro or pdb file that overrides the value for
  *bilayer* in the configuration file.
"""),
    ]



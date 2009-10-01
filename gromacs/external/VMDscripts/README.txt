.. $Id$

Using the VMD scripts
=====================

The VMD tcl scripts are currently not being used directly in code and
are purely experimental.

In order to use them in VMd (e.g. using the vmd remote module) one
will have to add the directory to the ``auto_path``. For instance,::

  # add local (autoloaded) scripts to the search path
  set auto_path [concat $env(HOME)/Library/python/GromacsWrapper/gromacs/external/VMDscripts $auto_path]

This is more complicated when installed as a zipped egg file; you will
have to extract the files manually and install elsewhere at the
moment. (Or do a developer install for the egg.)


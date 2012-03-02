.. -*- mode: rst -*-

=========
 INSTALL
=========

This document should help you to install the **GromacsWrapper**
package. The installation uses `setuptools`_ (also known as
``easy_install`` or "egg install"); if this is not available on your
system you can either let the installer download it automatically from
the internet (so just go to `Quick installation instructions`_) or
install it using your package manager, eg::

  aptitude install python-setuptools

or similar.

Please do not hesitate to contact `Oliver Beckstein`_ if problems
occur or if you have suggestions on how to improve the package or
these instructions.

.. _`Oliver Beckstein`: orbeckst@gmail.com
.. _setuptools: http://peak.telecommunity.com/DevCenter/setuptools


Quick installation instructions
===============================

If you have ``easy_install`` on your system you can directly install
from the interweb::

  easy_install -f https://github.com/orbeckst/GromacsWrapper/tags GromacsWrapper

This will automatically download and install the latest version.

Manual Download
===============

If your prefer to download manually, get the latest stable release
from

  https://github.com/orbeckst/GromacsWrapper/tags

and use any of the following methods (in increasing order of
complexity):

- From an egg install file, eg GromacsWrapper-0.1-py2.5.egg::

   easy_install GromacsWrapper-0.1-py2.5.egg

- From a tar ball, eg GromacsWrapper-0.1.tar.gz::

   easy_install GromacsWrapper-0.1.tar.gz

- From the unpacked source::

   tar -zxvf GromacsWrapper-0.1.tar.gz
   cd GromacsWrapper-0.1
   python setup.py install

See the `easy_install instructions`__ for explanation of the options
that allow you to install into non-standard places.

.. __: http://peak.telecommunity.com/DevCenter/EasyInstall#custom-installation-locations


Source code access
==================

The tar archive from https://github.com/orbeckst/GromacsWrapper/tags
contains a full source code distribution.

In order to follow code development you can also browse the code
**git** repository at http://github.com/orbeckst/GromacsWrapper or
clone the git repository from

   git://github.com/orbeckst/GromacsWrapper.git
 


Requirements
============

Python_ and Gromacs_ must be installed. ipython_ is very much
recommended. These packages might already be available through your local
package manager such as ``aptitude/apt``, ``yum``, ``yast``, ``fink`` or
``macports``. 

.. _Python: http://www.python.org
.. _Gromacs: http://www.gromacs.org
.. _ipython: http://ipython.scipy.org


System requirements
-------------------

Tested with python 2.5, 2.6 on Linux and Mac OS X. Earlier python
versions will likely fail.


Required python modules
-----------------------

The basic package makes use of numpy_ and can use matplotlib_ (in the
form of the ``pylab`` package). Only numpy_ is immediately required (and
automatically installed with ``easy_install``).

.. _numpy: http://numpy.scipy.org
.. _matplotlib: http://matplotlib.sourceforge.net/

For the :mod:`gromacs.analysis` library additional packages are required:

  =============  ==========  ==================================================
  package        version     source
  =============  ==========  ==================================================
  matplotlib     >=0.91.3    http://matplotlib.sourceforge.net/
  scipy                      http://www.scipy.org/
  RecSQL         >=0.3       http://sbcb.bioch.ox.ac.uk/oliver/software/RecSQL/
  =============  ==========  ==================================================

See `Installing all packages and requirements`_ for hints on how to
install these package.


Additional instructions
=======================

Installing *all* packages and requirements
------------------------------------------

If you want to make sure that ``easy_install`` also installs
requirements for optional modules then you will have to add the
additional requirement ``[analysis]`` to the command line. For a web
install this would look like ::

  easy_install -f https://github.com/orbeckst/GromacsWrapper/tags GromacsWrapper[analysis]

For installation from a downloaded source distribution ::

  easy_install GromacsWrapper-0.1.tar.gz[analysis]

or from within the unpacked source ::

  cd GromacsWrapper-0.1
  easy_install . GromacsWrapper[analysis]

In each case this will try to download additional packages for the
extra *analysis* module.

A common problem appears to be the error `Could not find matplotlib`_
as discussed below.


Troubleshooting
===============

For problems with ``easy_install`` please read the `User setuptools
instructions`__. 

.. __: http://peak.telecommunity.com/DevCenter/setuptools#what-your-users-should-know


Installing in non-standard locations
------------------------------------

Inform yourself about how to use ``easy_install`` to install packages
in `Custom Installation Locations`_.

For code hacking and development a *developer installation* is often
useful. In the unpacked source::

  python setup.py develop --install-dir python-lib-dir

where ``python-lib-dir`` must be on the ``PYTHONPATH``.

.. _Custom Installation Locations:
   http://peak.telecommunity.com/DevCenter/EasyInstall#custom-installation-locations


easy_install import error
-------------------------

Online installation can run into issues where the installation dies
with the error::

 ImportError: No module named ez_setup

If `EasyInstall Troubleshooting`_ does not help then try downloading
the source distribution package manually, unpack, and install from
inside with something like::

 python setup.py install

If this is still not working contact the author and complain.

.. _EasyInstall Troubleshooting:
   http://peak.telecommunity.com/DevCenter/EasyInstall#troubleshooting



Could not find matplotlib
-------------------------

Automatic downloading of ``matplotlib`` often fails::

   Searching for matplotlib>=0.91.3
   Reading http://pypi.python.org/simple/matplotlib/
   Reading http://matplotlib.sourceforge.net
   Reading https://sourceforge.net/project/showfiles.php?group_id=80706&package_id=278194
   Reading https://sourceforge.net/project/showfiles.php?group_id=80706&package_id=82474
   Reading http://sourceforge.net/project/showfiles.php?group_id=80706
   No local packages or download links found for matplotlib>=0.91.3
   error: Could not find suitable distribution for Requirement.parse('matplotlib>=0.91.3')

If automatic downloading of ``matplotlib`` fails then the best
approach is to install it through your package management
system. Search for "matplotlib" or "pylab" in the list of available packages.

If this is not an option then `download matplotlib`_ manually and
`install matplotlib manually`_ first. For example, ::

   wget http://kent.dl.sourceforge.net/sourceforge/matplotlib/matplotlib-0.98.5.3-py2.5-macosx-10.3-fat.egg \
        -O matplotlib-0.98.5.3-py2.5.egg

   easy_install matplotlib-0.98.5.3-py2.5.egg

Note that you should look at the `download matplotlib`_ page to get
the latest distribution. As highlighted in the `matplotlib
installation FAQ`_ it is important to rename the ``egg`` file (as done
in the example above).

Possibly the following installation from the source distribution
works, too::

   easy_install matplotlib-0.98.5.3.tar.gz

Once this has been accomplished, try the above installation
instructions again; ``easy_install`` should now pick up the newly
installed matplotlib.

.. _`download matplotlib`: 
     http://sourceforge.net/project/showfiles.php?group_id=80706
.. _`install matplotlib manually`: 
     http://matplotlib.sourceforge.net/users/installing.html
.. _`matplotlib installation FAQ`:
     http://matplotlib.sourceforge.net/faq/installing_faq.html#easy-install-from-egg

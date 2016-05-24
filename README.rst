.. -*- mode: rst -*-
.. The whole GromacsWrapper package is Copyright (c) 2009,2010,2011,2012 Oliver Beckstein,
.. except where noted otherwise.


========================
 README: GromacsWrapper
========================

|zenodo| |docs|

A primitive wrapper around the Gromacs tools until we have proper
python bindings. It also provides a small library (cook book) of
often-used recipes and an optional analysis module with plugins for
more complicated analysis tasks.

`Documentation`_ is mostly provided through the python doc strings and
available at http://gromacswrapper.readthedocs.org for recent releases.

The source code is available in the `GromacsWrapper git repository`_.

Please be aware that this is **alpha** software that most definitely
contains bugs. It is *your* responsibility to ensure that you are
running simulations with sensible parameters.

.. _Documentation: 
   http://gromacswrapper.readthedocs.org/en/latest/
.. _GromacsWrapper git repository:
   http://github.com/Becksteinlab/GromacsWrapper
.. |zenodo| image:: https://zenodo.org/badge/13219/Becksteinlab/GromacsWrapper.svg
   :target: https://zenodo.org/badge/latestdoi/13219/Becksteinlab/GromacsWrapper
   :alt: Latest release on zenodo (with DOI)
.. |docs| image:: https://readthedocs.org/projects/gromacswrapper/badge/?version=latest
   :target: http://gromacswrapper.readthedocs.org/en/latest/?badge=latest
   :alt: Documentation

Licence
=======

The **GromacsWrapper** package is made available under the terms of
the `GNU Public License v3`_ (or any higher version at your choice)
except as noted below. See the file COPYING for the licensing terms
for all modules.

The **vmd** module is made available under the `LGPL v3`_ (see COPYING
and COPYING.LESSER). **numkit** is provided under the "`Modified BSD
Licence`_" (as it contains some code from scipy_).

.. _GNU Public License v3: http://www.gnu.org/licenses/gpl.html
.. _LGPL v3: http://www.gnu.org/licenses/lgpl.html
.. _Modified BSD Licence: http://www.opensource.org/licenses/bsd-license.php
.. _scipy: http://www.scipy.org

The distribution contains third party software that is copyrighted by
the authors but distributed under licences compatible with this
package license. Where permitted and necessary, software/files were
modified to integrate with GromacsWrapper.


Installation
============

Releases
--------

The `latest version of GromacsWrapper from PyPi`_ and can be installed
with ::

  pip install GromacsWrapper

.. _`latest version of GromacsWrapper from PyPi`:
   https://pypi.python.org/pypi/GromacsWrapper

Development version
-------------------

The *develop* branch in the GitHub source repository generally
contains useful code but nevertheless, things can break in weird and
wonderful ways. Please report issues through the `Issue Tracker`_ and
mention that you used the *develop branch*.

To use the *development code base*:  checkout the ``develop`` branch::

   git clone https://github.com/Becksteinlab/GromacsWrapper.git
   cd GromacsWrapper
   git checkout -b develop origin/develop

and install ::

   python setup.py install




Download and Availability
=========================

The GromacsWrapper home page is
http://github.com/Becksteinlab/GromacsWrapper.  The latest release of the
package is being made available from https://github.com/Becksteinlab/GromacsWrapper/releases

You can also clone the `GromacsWrapper git repository`_ or fork for
your own development::

  git clone git://github.com/Becksteinlab/GromacsWrapper.git


Reporting Bugs and Contributing to GromacsWrapper
=================================================

Please use the `Issue Tracker`_ to report bugs, installation problems,
and feature requests.

**Pull requests** for bug fixes and enhancements are very welcome.

.. _Issue Tracker: http://github.com/Becksteinlab/GromacsWrapper/issues

Building Documentation
======================

Install Sphinx::

   apt-get install python-sphinx

and compile::

   cd package/doc/sphinx
   make html



Citing
======

|zenodo|

GromacsWrapper was written by Oliver Beckstein with contributions from
many other people. Please see the file AUTHORS_ for all the names.

If you find this package useful and use it in published work I'd be
grateful if it was acknowledged in text as

  "... used GromacsWrapper (Oliver Beckstein et al,
  http://github.com/Becksteinlab/GromacsWrapper doi: 10.5281/zenodo.17901)"

or in the Acknowledgements section.

Thank you.

.. _AUTHORS:
   https://raw.githubusercontent.com/Becksteinlab/GromacsWrapper/develop/AUTHORS


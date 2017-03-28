.. -*- mode: rst -*-
.. The whole GromacsWrapper package is Copyright (c) 2009,2010,2011,2012 Oliver Beckstein,
.. except where noted otherwise.


========================
 README: GromacsWrapper
========================

|zenodo|

A primitive wrapper around the Gromacs tools until we have proper
python bindings. It also provides a small library (cook book) of
often-used recipes and an optional analysis module with plugins for
more complicated analysis tasks.


See :doc:`INSTALL` for installation instructions. `Documentation`_ is
mostly provided through the python doc strings. See `Download and
Availability`_ for download instructions if the instructions in
:doc:`INSTALL` are not sufficient.

The source code is also available in the `GromacsWrapper git
repository`_.

Please be aware that this is **alpha** software that most definitely
contains bugs. It is *your* responsibility to ensure that you are
running simulations with sensible parameters.


.. _Documentation: 
   http://gromacswrapper.readthedocs.org
.. _GromacsWrapper git repository:
   http://github.com/Becksteinlab/GromacsWrapper
.. |zenodo| image:: https://zenodo.org/badge/13219/Becksteinlab/GromacsWrapper.svg
   :target: https://zenodo.org/badge/latestdoi/13219/Becksteinlab/GromacsWrapper
   :alt: Latest release on zenodo (with DOI)

Licence
=======

The **GromacsWrapper** package is made available under the terms of
the `GNU Public License v3`_ (or any higher version at your choice)
except as noted below. See the file COPYING for the licensing terms
for all modules.

.. _GNU Public License v3: http://www.gnu.org/licenses/gpl.html


Included Software
=================

The distribution contains third party software that is copyrighted by
the authors but distributed under licences compatible with this
package license. Where permitted and necessary, software/files were
modified to integrate with GromacsWrapper.

In case of problems please direct error reports through the
GromacsWrapper `Issue Tracker`_ as these bugs might not have been
present in the original software or files.

Included third party content:

``gromacs.fileformats.preprocess``
  - The preprocessor is based on ``pypreprocessor.py`` from
    `pypreprocessor`_, release 0.4.0.
  - copyright: (c) 2010 Evan Plaice
  - license: MIT
  - https://github.com/evanplaice/pypreprocessor

  .. _pypreprocessor: https://github.com/evanplaice/pypreprocessor


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



Download and Availability
=========================

The GromacsWrapper home page is
http://github.com/Becksteinlab/GromacsWrapper.  The latest release of the
package is being made available from https://github.com/Becksteinlab/GromacsWrapper/releases

You can also clone the `GromacsWrapper git repository`_ or fork for
your own development::

  git clone git://github.com/Becksteinlab/GromacsWrapper.git



Contact
=======

Please use the `Issue Tracker`_ to report bugs, installation problems,
and feature requests (mention ``@orbeckst`` in the issue report).

.. _Issue Tracker: http://github.com/Becksteinlab/GromacsWrapper/issues

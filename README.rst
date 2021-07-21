.. -*- mode: rst, coding: utf-8 -*-
.. The whole GromacsWrapper package is Copyright (c) 2009-2018 Oliver
.. Beckstein and AUTHORS except where noted otherwise.


========================
 README: GromacsWrapper
========================

|build| |cov| |docs| |zenodo| |PRsWelcome|

A primitive Python wrapper around the Gromacs_ tools. The library is
tested with Gromacs 4.6.5, 2018.x, 2019.x (and 2016.x also works) and
is tested with Python 2.7 and 3.6--3.9.

GromacsWrapper also provides a small library (cook book) of often-used
recipes and helper functions to set up MD simulations.

`Documentation`_ is mostly provided through the python doc strings and
available at https://gromacswrapper.readthedocs.org for recent releases.

The source code is available in the `GromacsWrapper git repository`_.

Please be aware that this is **beta** software that most definitely
contains bugs. It is *your* responsibility to ensure that you are
running simulations with sensible parameters.

.. _Gromacs: http://www.gromacs.org
.. _Documentation: 
   https://gromacswrapper.readthedocs.org/en/latest/
.. _GromacsWrapper git repository:
   https://github.com/Becksteinlab/GromacsWrapper
.. |build| image:: https://github.com/Becksteinlab/GromacsWrapper/actions/workflows/ci.yaml/badge.svg?branch=master
   :target: https://github.com/Becksteinlab/GromacsWrapper/actions/workflows/ci.yaml
   :alt: Build Status	     
.. |cov| image:: https://codecov.io/gh/Becksteinlab/GromacsWrapper/badge.svg
   :target: https://codecov.io/gh/Becksteinlab/GromacsWrapper
   :alt: Code Coverage
   :scale: 100%
.. |zenodo| image:: https://zenodo.org/badge/13219/Becksteinlab/GromacsWrapper.svg
   :target: https://zenodo.org/badge/latestdoi/13219/Becksteinlab/GromacsWrapper
   :alt: Latest release on zenodo (with DOI)
.. |docs| image:: https://readthedocs.org/projects/gromacswrapper/badge/?version=latest
   :target: https://gromacswrapper.readthedocs.org/en/latest/?badge=latest
   :alt: Documentation
.. |PRsWelcome| image:: https://img.shields.io/badge/PRs-welcome-brightgreen.svg
   :target: http://makeapullrequest.com
   :alt: PRs Welcome!

Quick Start
===========

Given a PDB file ``1ake.pdb``, set up and run a simple simulation (assuming
you have all other input files at hand such as the MDP files)::

  >>> import gromacs
  >>> print(gromacs.release)
  2018.2
  >>> help(gromacs.pdb2gmx)
  DESCRIPTION

  gmx pdb2gmx reads a .pdb (or .gro) file, reads some database files,
  adds hydrogens to the molecules and generates coordinates in GROMACS
  ...
  ...
  OPTIONS

  Options to specify input files:

  -f      [<.gro/.g96/...>]  (eiwit.pdb)
            Structure file: gro g96 pdb brk ent esp tpr
  ...
  ...
  >>> gromacs.pdb2gmx(f="1ake.pdb", o="protein.gro", p="topol.top",
  ...                 ff="oplsaa", water="tip4p")
  >>> gromacs.editconf(f="protein.gro", o="boxed.gro",
  ...                  bt="dodecahedron", d=1.5, princ=True,
  ...                  input="Protein")
  >>> gromacs.solvate(cp="boxed.gro", cs="tip4p", p="topol.top",
  ...                 o="solvated.gro")
  >>> gromacs.grompp(f="emin.mdp", c="solvated.gro", p="topol.top",
  ...                o="emin.tpr")
  >>> gromacs.mdrun(v=True, deffnm="emin")
  >>> gromacs.grompp(f="md.mdp", c="emin.gro", p="topol.top", o="md.tpr")
  >>> gromacs.mdrun(v=True, deffnm="md")


	 
License
=======

The **GromacsWrapper** package is made available under the terms of
the `GNU Public License v3`_ (or any higher version at your choice)
except as noted below. See the file COPYING for the licensing terms
for all modules.

.. _GNU Public License v3: http://www.gnu.org/licenses/gpl.html

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
   https://pypi.org/project/GromacsWrapper/

Development version
-------------------

The *master* branch in the GitHub source repository generally
contains useful code but nevertheless, things can break in weird and
wonderful ways. Please report issues through the `Issue Tracker`_ and
mention that you used the *develop branch*.

To use the *development code base*:  checkout the ``master`` branch::

   git clone https://github.com/Becksteinlab/GromacsWrapper.git
   cd GromacsWrapper

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

Questions
=========

Please ask questions in the `Issue Tracker`_ (instead of private email).


Reporting Bugs and Contributing to GromacsWrapper
=================================================

Please use the `Issue Tracker`_ to report bugs, installation problems,
and feature requests.

**Pull requests** for bug fixes and enhancements are very welcome. See http://makeapullrequest.com for a 
general introduction on how make a pull request and contribute to open source projects.

.. _Issue Tracker: http://github.com/Becksteinlab/GromacsWrapper/issues


Building Documentation
======================

Install Sphinx::

   pip install sphinx

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
  https://github.com/Becksteinlab/GromacsWrapper doi: 10.5281/zenodo.17901)"

or in the Acknowledgements section.

Thank you.

.. _AUTHORS:
   https://raw.githubusercontent.com/Becksteinlab/GromacsWrapper/master/AUTHORS


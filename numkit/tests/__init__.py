# -*- coding: utf-8 -*-
# numkit test cases
# Part of GromacsWrapper
# Copyright (c) Oliver Beckstein <orbeckst@gmail.com>
# Published under the Modified BSD Licence.

"""
=========================
 Test cases for numkit
=========================

We are using the NumPy_ testing frame work; thus, numpy *must* be
installed for the tests to run at all.

Run all the tests with

   >>> import numkit.tests
   >>> numkit.tests.test(label='full')

Additional information is displayed at a higher verbosity level (the default is 1):

   >>> numkit.tests.test(label='fast', verbose=3)

Note that if no tests are being run then one might have to run the
tests with the ``--exe`` flag

   >>> numkit.tests.test(label='fast', extra_argv=['--exe'])

(This happens when python files are installed with the executable bit set. By
default the nose_ testing framework refuses to use those files and must be
encouraged to do so with the ``--exe`` switch.)

See `nose commandline options`_ for additional options that can be used; for
instance, code coverage can also be checked:

  >>> numkit.tests.test(label='full', extra_argv=['--exe', '--with-coverage'])


Writing test cases
==================

The unittests use the :mod:`unittest` module together with nose_. See the
examples in the ``numkit/tests`` directory.

The `SciPy testing guidelines`_ are a good howto for writing test cases,
especially as we are directly using this framework (imported from numpy).


.. _nose: 
   http://somethingaboutorange.com/mrl/projects/nose/0.11.3/index.html
.. _nose commandline options:
   http://somethingaboutorange.com/mrl/projects/nose/0.11.3/usage.html#extended-usage
.. _SciPy testing guidelines: 
   http://projects.scipy.org/numpy/wiki/TestingGuidelines#id11
"""

try:
    from numpy.testing import Tester
    test = Tester().test
except ImportError:
    raise ImportError("""numpy>=1.3  is required to run the test suite. Please install it first. """
                      """(For example, try "easy_install 'numpy>=1.3'").""")


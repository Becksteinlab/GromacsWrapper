# GromacsWrapper
# Copyright (c) 2009-2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
:mod:`gromacs` -- GromacsWrapper Package Overview
=================================================

**GromacsWrapper** (package :mod:`gromacs`) is a thin shell around the `Gromacs`_
tools for light-weight integration into python scripts or interactive use in
`ipython`_.

.. _`Gromacs`: http://www.gromacs.org
.. _`ipython`: http://ipython.scipy.org


Modules
-------

:mod:`gromacs`
     The top level module contains all gromacs tools; each tool can be
     run directly or queried for its documentation. It also defines
     the root logger class (name *gromacs* by default).

:mod:`gromacs.config`
     Configuration options. Not really used much at the moment.

:mod:`gromacs.cbook`
     The Gromacs cook book contains typical applications of the tools. In many
     cases this not more than just an often-used combination of parameters for
     a tool.

:mod:`gromacs.tools`
     Contains classes that wrap the gromacs tools. They are automatically
     generated from the list of tools in :data:`gromacs.tools.gmx_tools`.

:mod:`gromacs.formats`
     Classes to represent data files in various formats such as
     xmgrace graphs. The classes allow reading and writing and for
     graphs, also plotting of the data.

:mod:`gromacs.utilities`
     Convenience functions and mixin-classes that are used as helpers
     in other modules.

:mod:`gromacs.setup`
     Functions to set up a MD simulation, containing tasks such as solvation
     and adding ions, energy minimizqtion, MD with position-restraints, and
     equilibrium MD.

:mod:`gromacs.qsub`
     Functions to handle batch submission queuing systems.

:mod:`gromacs.run`
     Classes to run :program:`mdrun` in various way, including on
     multiprocessor systems.

:mod:`gromacs.analysis`
     A package that collects whole analysis tasks. It uses the
     :mod:`gromacs` package but is otherwise only loosely coupled with
     the rest. At the moment it only contains the infrastructure and
     an example application. See the package documentation.


Examples
--------

The following examples should simply convey the flavour of using the
package. See the individual modules for more examples.

Getting help
............

In python::

   help(gromacs.g_dist)
   gromacs.g_dist.help()
   gromacs.g_dist.help(long=True)

In ``ipython``::

   gromacs.g_dist ?


Simple usage
............

Gromacs flags are given as python keyword arguments::

   gromacs.g_dist(v=True, s='topol.tpr', f='md.xtc', o='dist.xvg', dist=1.2)

Input to stdin of the command can be supplied::

   gromacs.make_ndx(f='topol.tpr', o='md.ndx',
                    input=('keep "SOL"', '"SOL" | r NA | r CL', 'name 2 solvent', 'q'))

Output of the command can be caught in a variable and analyzed::

   rc, output, junk = gromacs.grompp(..., stdout=False)        # collects command output
   for line in output.split('\\n'):
       line = line.strip()
       if line.startswith('System has non-zero total charge:'):
             qtot = float(line[34:])
             break

(See :func:`gromacs.cbook.grompp_qtot` for a more robust implementation of this
application.)


Warnings and Exceptions
-----------------------

A number of package-specific exceptions (:exc:`GromacsError`) and
warnings (:exc:`GromacsFailureWarning`, :exc:`GromacsImportWarning`, 
:exc:`GromacsValueWarning`, :exc:`AutoCorrectionWarning`,
:exc:`BadParameterWarning`) can be raised.

If you want to stop execution at, for instance, a :exc:`AutoCorrectionWarning` or
:exc:`BadParameterWarning` then use the python :mod:`warnings` filter::

  import warnings
  warnings.simplefilter('error', gromacs.AutoCorrectionWarning)
  warnings.simplefilter('error', gromacs.BadParameterWarning)

This will make python raise an exception instead of moving on. The default is
to always report, eg::

  warnings.simplefilter('always', gromacs.BadParameterWarning)

The following *exceptions* are defined:

.. autoexception:: GromacsError
.. autoexception:: MissingDataError
.. autoexception:: ParseError

The following *warnings* are defined:

.. autoexception:: GromacsFailureWarning
.. autoexception:: GromacsImportWarning
.. autoexception:: GromacsValueWarning
.. autoexception:: AutoCorrectionWarning
.. autoexception:: BadParameterWarning
.. autoexception:: MissingDataWarning
.. autoexception:: UsageWarning
.. autoexception:: LowAccuracyWarning


Logging
-------

The library uses python's logging_ module to keep a history of what it has been
doing. In particular, every wrapped Gromacs command logs its command line
(including piped input) to the log file (configured in
:data:`gromacs.config.logfilename`). This facilitates debugging or simple
re-use of command lines for very quick and dirty work. The logging facilty
appends to the log file and time-stamps every entry. See :mod:`gromacs.config`
for more details on configuration.

It is also possible to capture output from Gromacs commands in a file
instead of displaying it on screen, as described under
:ref:`input-output-label`.

.. _logging: http://docs.python.org/library/logging.html

Version
-------

The package version can be queried with the :func:`gromacs.get_version` function.

.. autofunction:: get_version
.. autofunction:: get_version_tuple

If the package was installed from a development version, the patch
level will have the string "-dev" affixed to distinguish it from a
release.
"""
from __future__ import absolute_import
__docformat__ = "restructuredtext en"

from .version import VERSION, RELEASE, get_version, get_version_tuple

# __all__ is extended with all gromacs command instances later
__all__ = ['config', 'tools', 'cbook', 'fileformats']

from . import fileformats

# Note: analysis not imported by default (requires additional pre-requisites)

import warnings
from .exceptions import (GromacsError, MissingDataError, ParseError,
                         GromacsFailureWarning, GromacsImportWarning,
                         GromacsValueWarning, AutoCorrectionWarning,
                         BadParameterWarning, MissingDataWarning,
                         UsageWarning, LowAccuracyWarning)


# Import configuration before anything else
from . import config


import logging
# NOTE: logging is still iffy; when I reload I add a new logger each
# time and output is repeated for each reload. Probably should heed
# the advice on logging and libraries in
# http://docs.python.org/library/logging.html?#configuring-logging-for-a-library
class NullHandler(logging.Handler):
    def emit(self, record):
        pass

# default silent logger --- just here for illustration; below we
# we get a proper logger from log.create()
h = NullHandler()
logging.getLogger("gromacs").addHandler(h)
del h

# The top level logger of the library is named 'gromacs'.
# Note that we are configuring this logger with console output. If the root logger also
# does this then we will get two output lines to the console. We'll live with this because
# this is a simple convenience library and most people will not bother
# with a logger (I think...)
#
# In modules that use loggers get a logger like so:
#     import logging
#     logger = logging.getLogger('gromacs.MODULENAME')

def start_logging(logfile="gromacs.log"):
    """Start logging of messages to file and console.

    The default logfile is named ``gromacs.log`` and messages are
    logged with the tag *gromacs*.
    """
    import log
    log.create("gromacs", logfile=logfile)
    logging.getLogger("gromacs").info("GromacsWrapper %s STARTED logging to %r", get_version(), logfile)

def stop_logging():
    """Stop logging to logfile and console."""
    import log
    logger = logging.getLogger("gromacs")
    logger.info("GromacsWrapper %s STOPPED logging", get_version())
    log.clear_handlers(logger)  # this _should_ do the job...


# Add gromacs command **instances** to the top level.
# These serve as the equivalence of running commands in the shell.
# (Note that each gromacs command is actually run when the instance is
# created in order to gather the documentation string.)
from . import tools

# Ignore warnings from a few programs that do not produce
# documentation when run with '-h' (only applies when the default for
# failuremode of core.GromacsCommand is changed to 'warn')
warnings.simplefilter("ignore", GromacsFailureWarning)
_have_g_commands = []
_missing_g_commands = []
for clsname, cls in tools.registry.items():
    name = clsname[0].lower() + clsname[1:]    # instances should start with lower case
    try:
        globals()[name] = cls()                # add instance of command for immediate use
        _have_g_commands.append(name)
    except:
        _missing_g_commands.append(name)
warnings.simplefilter("always", GromacsFailureWarning)
warnings.simplefilter("always", GromacsImportWarning)

_have_g_commands.sort()
_missing_g_commands.sort()
if len(_missing_g_commands) > 0:
    warnings.warn("Some Gromacs commands were NOT found; "
                  "maybe source GMXRC first? The following are missing:\n%r\n" % _missing_g_commands,
                  category=GromacsImportWarning)

del name, cls, clsname

# get ALL active command instances with 'from gromacs import *'
__all__.extend(_have_g_commands)


# cbook should come after the whole of init as it relies on command
# instances in the top level name space
try:
    from . import cbook
except OSError, err:
    warnings.warn("Some Gromacs commands were NOT found when importing gromacs.cbook:\n"+str(err),
                  category=GromacsImportWarning)

# convenience functions for warnings

less_important_warnings = ['AutoCorrectionWarning', 'UsageWarning']

def filter_gromacs_warnings(action, categories=None):
    """Set the :meth:`warnings.simplefilter` to *action*.

    *categories* must be a list of warning classes or strings.
    ``None`` selects the defaults,  :data:`gromacs.less_important_warnings`.
    """

    if categories is None:
        categories = less_important_warnings
    for c in categories:
        try:
            w = globals()[c]
        except KeyError:
            w = c
        if not issubclass(w, Warning):
            raise TypeError("%r is neither a Warning nor the name of a Gromacs warning." % c)
        warnings.simplefilter(action, category=w)

def disable_gromacs_warnings(categories=None):
    """Disable ("ignore") specified warnings from the gromacs package.

    *categories* must be a list of warning classes or strings.
    ``None`` selects the defaults.

    """
    filter_gromacs_warnings('ignore', categories=categories)

def enable_gromacs_warnings(categories=None):
    """Enable ("always") specified warnings from the gromacs package.

    *categories* must be a list of warning classes or strings.
    ``None`` selects the defaults, :data:`gromacs._less_important_warnings`.

    """
    filter_gromacs_warnings('always', categories=categories)


# define the testing framework
from numpy.testing.nosetester import NoseTester

test = NoseTester().test
del NoseTester

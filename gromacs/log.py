# log.py
# logging for Gromacs module
# Copyright (c) 2010 Oliver Beckstein
# Published under the GNU Public Licence v3

"""
:mod:`gromacs.log` -- setting up logging
========================================

Configure logging for GromacsWrapper. Import this module if logging is
desired in application code.

Logging to a file and the console is set up by default as described
under `logging to multiple destinations`_.

The top level logger of the library is named *gromacs* by convention;
a simple logger that writes to the console and logfile can be created
with the :func:`log.create` function. This only has to be done
*once*. For convenience, the default GromacsWrapper logger can be
created with :func:`gromacs.start_logger`::

 import gromacs
 gromacs.start_logger()

Once this has been done, GromacsWrapper will write messages to the
logfile (named :file:`gromacs.log` by default although this can be
changed with the optional argument to :func:`~gromacs.start_logger`).

Any code can log to the gromacs logger by using ::

 import logging
 logger = logging.getLogger('gromacs.MODULENAME')

 # use the logger, for example at info level:
 logger.info("Starting task ...")

The important point is that the name of the logger begins with
"gromacs.".

.. _logging to multiple destinations:
   http://docs.python.org/library/logging.html?#logging-to-multiple-destinations

.. SeeAlso:: The :mod:`logging` module in the standard library contains
             in depth documentation about using logging.


Convenience functions
---------------------

Two convenience functions at the top level make it easy to start and
stop the default *gromacs* logger.

.. autofunction:: gromacs.start_logging
.. autofunction:: gromacs.stop_logging


Other functions and classes for logging purposes
------------------------------------------------

.. autogenerated, see Online Docs

"""

import logging


def create(logger_name, logfile="gromacs.log"):
    """Create a top level logger.

    - The file logger logs everything (including DEBUG).
    - The console logger only logs INFO and above.

    Logging to a file and the console.

    See http://docs.python.org/library/logging.html?#logging-to-multiple-destinations

    The top level logger of the library is named 'gromacs'.  Note that
    we are configuring this logger with console output. If the root
    logger also does this then we will get two output lines to the
    console. We'll live with this because this is a simple
    convenience library...
    """

    logger = logging.getLogger(logger_name)

    logger.setLevel(logging.DEBUG)

    logfile = logging.FileHandler(logfile)
    logfile_formatter = logging.Formatter(
        "%(asctime)s %(name)-12s %(levelname)-8s %(message)s"
    )
    logfile.setFormatter(logfile_formatter)
    logger.addHandler(logfile)

    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter("%(name)-12s: %(levelname)-8s %(message)s")
    console.setFormatter(formatter)

    logger.addHandler(console)

    return logger


def clear_handlers(logger):
    """clean out handlers in the library top level logger

    (only important for reload/debug cycles...)
    """
    for h in logger.handlers:
        logger.removeHandler(h)


class NullHandler(logging.Handler):
    """Silent Handler.

    Useful as a default::
      h = NullHandler()
      logging.getLogger("gromacs").addHandler(h)
      del h
    """

    def emit(self, record):
        pass

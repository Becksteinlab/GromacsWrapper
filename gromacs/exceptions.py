# GromacsWrapper
# Copyright (c) 2009-2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

# exceptions and warnings

class GromacsError(EnvironmentError):
    """Error raised when a gromacs tool fails.

    Returns error code in the errno attribute and a string in strerror.
    # TODO: return status code and possibly error message
    """

class MissingDataError(Exception):
    """Error raised when prerequisite data are not available.

    For analysis with :class:`gromacs.analysis.core.Simulation` this typically
    means that the :meth:`~gromacs.analysis.core.Simulation.analyze` method has
    to be run first.
    """

class ParseError(Exception):
    """Error raised when parsing of a file failed."""

class GromacsFailureWarning(Warning):
    """Warning about failure of a Gromacs tool."""

class GromacsImportWarning(ImportWarning):
    """Warns about problems with using a gromacs tool."""

class GromacsValueWarning(Warning):
    """Warns about problems with the value of an option or variable."""

class AutoCorrectionWarning(Warning):
    """Warns about cases when the code is choosing new values automatically."""

class BadParameterWarning(Warning):
    """Warns if some parameters or variables are unlikely to be appropriate or correct."""

class MissingDataWarning(Warning):
    """Warns when prerequisite data/files are not available."""

class UsageWarning(Warning):
    """Warns if usage is unexpected/documentation ambiguous."""

class LowAccuracyWarning(Warning):
    """Warns that results may possibly have low accuracy."""

import warnings
# These warnings should always be displayed because other parameters
# can have changed, eg during interactive use.
for w in (AutoCorrectionWarning, BadParameterWarning, UsageWarning,
          GromacsFailureWarning, GromacsValueWarning, LowAccuracyWarning):
    warnings.simplefilter('always', category=w)
del w

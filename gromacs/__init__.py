# $Id$
"""\
Gromacs shell
=============

A thin shell around the Gromacs tools for light-weight integration into python
scripts or interactive use in ``ipython``.


"""

__all__ = ['tools', 'analysis']

class GromacsError(EnvironmentError):
    """Error raised when a gromacs tool fails.

    Returns error code in the errno attribute and a string in strerror.
    """
    # TODO: return status code and possibly error message
    pass

class GromacsFailureWarning(Warning):
    """Warning about failure of a Gromacs tool."""
    pass

class GromacsImportWarning(ImportWarning):
    """Warns about problems with using a gromacs tool."""
    pass

# add gromacs command **instances** (note that each gromacs command is
# actually run when the instance is created in order to gather the
# documentation string.)
import tools

# Ignore warnings from a few programs that do not produce
# documentation when run with '-h' (only applies when the default for
# failuremode of core.GromacsCommand is changed to 'warn')
import warnings
warnings.simplefilter("ignore", GromacsFailureWarning)
_have_g_commands = []
_missing_g_commands = []
for clsname, cls in tools.registry.items():
    name = clsname[0].lower() + clsname[1:]    # instances should start lower case
    try:
        locals()[name] = cls()                 # add instance of command for immediate use
        _have_g_commands.append(name)
    except (GromacsError, OSError):            # ignore missing commands/missing -h for doc extraction
        _missing_g_commands.append(name)
warnings.simplefilter("default", GromacsFailureWarning)

warnings.simplefilter("default", GromacsImportWarning)
_have_g_commands.sort()
_missing_g_commands.sort()
if len(_missing_g_commands) > 0:
    warnings.warn("Some Gromacs commands were NOT found; "
                  "maybe source GMXRC first? The following are missing:\n%r\n" % _missing_g_commands,
                  category=GromacsImportWarning)

del name, cls, clsname

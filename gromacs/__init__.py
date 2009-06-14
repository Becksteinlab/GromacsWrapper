# $Id$
"""\
Gromacs shell
=============

A thin shell around the Gromacs tools for light-weight integration into python
scripts or interactive use in ``ipython``.


"""

__all__ = ['tools', 'analysis']

# add gromacs command **instances** (note that each gromacs command is
# actually run when the instance is created in order to gather the
# documentation string.)
import tools

# ignore warnings from a few programs that do not produce documentation when
# run with '-h'
import warnings
warnings.simplefilter("ignore", RuntimeWarning)
_have_g_commands = []
_missing_g_commands = []
for clsname, cls in tools.registry.items():
    name = clsname[0].lower() + clsname[1:]    # instances should start lower case
    try:
        locals()[name] = cls()                 # add instance of command for immediate use
        _have_g_commands.append(name)
    except OSError,err:
        _missing_g_commands.append(name)
warnings.simplefilter("default", RuntimeWarning)
_have_g_commands.sort()
_missing_g_commands.sort()
if len(_missing_g_commands) > 0:
    warnings.warn("Some Gromacs commands were NOT found; "
                  "maybe source GMXRC first? The following are missing:\n%r\n" % _missing_g_commands,
                  category=RuntimeWarning)

del name, cls, clsname

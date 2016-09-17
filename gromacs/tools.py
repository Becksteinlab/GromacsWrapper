# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.
"""
:mod:`gromacs.tools` -- Gromacs commands classes
================================================

A Gromacs command class produces an instance of a Gromacs tool command (
:class:`gromacs.core.GromacsCommand`), any argument or keyword argument
supplied will be used as default values for when the command is run.

Classes has the same name of the corresponding Gromacs tool with the first
letter capitalized and dot and dashes replaced by underscores to make it a
valid python identifier. Gromacs 5 tools are also aliased to their Gromacs 4
tool names (e.g, `sasa` to `g_sas`) for backwards compatibility.

The list of tools to be loaded is configured with the ``tools`` and ``groups``
options of the ``~/.gromacswrapper.cfg`` file. Guesses are made if these
options are not provided.

In the following example we create two instances of the
:class:`gromacs.tools.Trjconv` command (which runs the Gromacs ``trjconv``
command)::

  from gromacs.tools import Trjconv

  trjconv = tools.Trjconv()
  trjconv_compact = tools.Trjconv(ur='compact', center=True, boxcenter='tric', pbc='mol',
                                  input=('protein','system'))

The first one, ``trjconv``, behaves as the standard commandline tool but the
second one, ``trjconv_compact``, will by default create a compact
representation of the input data by taking into account the shape of the unit
cell. Of course, the same effect can be obtained by providing the corresponding
arguments to ``trjconv`` but by naming the more specific command differently
one can easily build up a library of small tools that will solve a specific,
repeatedly encountered problem reliably. This is particularly helpful when doing
interactive work.

Multi index
-----------

It is possible to extend the tool commands and patch in additional
functionality. For example, the :class:`GromacsCommandMultiIndex` class makes a
command accept multiple index files and concatenates them on the fly; the
behaviour mimics Gromacs' "multi-file" input that has not yet been enabled for
all tools.

.. autoclass:: GromacsCommandMultiIndex
.. autofunction:: merge_ndx

Helpers
-------

.. autofunction:: tool_factory
.. autofunction:: load_v4_tools
.. autofunction:: load_v5_tools
.. autofunction:: find_executables
.. autofunction:: make_valid_identifier
.. autoexception:: GromacsToolLoadingError

Gromacs tools
-------------

"""
from __future__ import absolute_import

import os.path
import tempfile
import subprocess
import atexit
import logging

from . import config
from .core import GromacsCommand

logger = logging.getLogger("gromacs.tools")

V4TOOLS = ("g_cluster", "g_dyndom", "g_mdmat", "g_principal", "g_select",
           "g_wham", "mdrun", "do_dssp", "g_clustsize", "g_enemat", "g_membed",
           "g_protonate", "g_sgangle", "g_wheel", "mdrun_d", "editconf",
           "g_confrms", "g_energy", "g_mindist", "g_rama", "g_sham", "g_x2top",
           "mk_angndx", "eneconv", "g_covar", "g_filter", "g_morph", "g_rdf",
           "g_sigeps", "genbox", "pdb2gmx", "g_anadock", "g_current",
           "g_gyrate", "g_msd", "g_sorient", "genconf", "g_anaeig", "g_density",
           "g_h2order", "g_nmeig", "g_rms", "g_spatial", "genion", "tpbconv",
           "g_analyze", "g_densmap", "g_hbond", "g_nmens", "g_rmsdist",
           "g_spol", "genrestr", "trjcat", "g_angle", "g_dielectric", "g_helix",
           "g_nmtraj", "g_rmsf", "g_tcaf", "gmxcheck", "trjconv", "g_bar",
           "g_dih", "g_helixorient", "g_order", "g_rotacf", "g_traj", "gmxdump",
           "trjorder", "g_bond", "g_dipoles", "g_kinetics", "g_pme_error",
           "g_rotmat", "g_tune_pme", "grompp", "g_bundle", "g_disre", "g_lie",
           "g_polystat", "g_saltbr", "g_vanhove", "make_edi", "xpm2ps", "g_chi",
           "g_dist", "g_luck", "g_potential", "g_sas", "g_velacc", "make_ndx")


#: dict of names in Gromacs 5 that correspond to an equivalent tool in
#: in Gromacs 4. The names are literal Gromacs names.
NAMES5TO4 = {
    # same name in both versions
    'grompp': 'grompp',
    'eneconv': 'eneconv',
    'editconf': 'editconf',
    'pdb2gmx': 'pdb2gmx',
    'trjcat': 'trjcat',
    'trjconv': 'trjconv',
    'trjorder': 'trjorder',
    'xpm2ps': 'xpm2ps',
    'mdrun': 'mdrun',
    'make_ndx': 'make_ndx',
    'make_edi': 'make_edi',
    'genrestr': 'genrestr',
    'genion': 'genion',
    'genconf': 'genconf',
    'do_dssp': 'do_dssp',

    # changed names
    'convert-tpr': 'tpbconv',
    'dump': 'gmxdump',
    'check': 'gmxcheck',
    'solvate': 'genbox',
    'distance': 'g_dist',
    'sasa': 'g_sas',
    'gangle': 'g_sgangle'
}


class GromacsToolLoadingError(Exception):
    """Raised when no Gromacs tool could be found."""


class GromacsCommandMultiIndex(GromacsCommand):
    """ Command class that accept multiple index files.

    It works combining multiple index files into a single temporary one so
    that tools that do not (yet) support multi index files as input can be
    used as if they did.

    It creates a new file only if multiple index files are supplied.
    """
    def __init__(self, **kwargs):
        kwargs = self._fake_multi_ndx(**kwargs)
        super(GromacsCommandMultiIndex, self).__init__(**kwargs)

    def run(self,*args,**kwargs):
        kwargs = self._fake_multi_ndx(**kwargs)
        return super(GromacsCommandMultiIndex, self).run(*args, **kwargs)

    def _fake_multi_ndx(self, **kwargs):
        ndx = kwargs.get('n')
        if not (ndx is None or type(ndx) is basestring) and \
           len(ndx) > 1 and 's' in kwargs:
            ndx.append(kwargs.get('s'))
            kwargs['n'] = merge_ndx(*ndx)
        return kwargs


def tool_factory(clsname, name, driver, base=GromacsCommand):
    """ Factory for GromacsCommand derived types. """
    clsdict = {
        'command_name': name,
        'driver': driver,
        '__doc__': property(base._get_gmx_docs)
    }
    return type(clsname, (base,), clsdict)


def make_valid_identifier(name):
    """ Turns tool names into valid identifiers.

    :param name: tool name
    :return: valid identifier
    """
    return name.replace('-', '_').capitalize()


def find_executables(path):
    """ Find executables in a path.

    Searches executables in a directory excluding some know commands
    unusable with GromacsWrapper.

    :param path: dirname to search for
    :return: list of executables
    """
    execs = []
    for exe in os.listdir(path):
        fullexe = os.path.join(path, exe)
        if (os.access(fullexe, os.X_OK) and not os.path.isdir(fullexe) and
             exe not in ['GMXRC', 'GMXRC.bash', 'GMXRC.csh', 'GMXRC.zsh',
                         'demux.pl', 'xplor2gmx.pl']):
            execs.append(exe)
    return execs


def load_v5_tools():
    """ Load Gromacs 5.x tools automatically using some heuristic.

    Tries to load tools (1) using the driver from configured groups (2) and
    falls back to automatic detection from ``GMXBIN`` (3) then to rough guesses.

    In all cases the command ``gmx help`` is ran to get all tools available.

    :return: dict mapping tool names to GromacsCommand classes
    """
    logger.debug("Loading v5 tools...")

    drivers = config.get_tool_names()

    if len(drivers) == 0 and 'GMXBIN' in os.environ:
        drivers = find_executables(os.environ['GMXBIN'])

    if len(drivers) == 0 or len(drivers) > 4:
        drivers = ['gmx', 'gmx_d', 'gmx_mpi', 'gmx_mpi_d']

    tools = {}
    for driver in drivers:
        try:
            out = subprocess.check_output([driver, '-quiet', 'help',
                                           'commands'])
            for line in str(out).encode('ascii').splitlines()[5:-1]:
                if line[4] != ' ':
                    name = line[4:line.index(' ', 4)]
                    fancy = make_valid_identifier(name)
                    suffix = driver.partition('_')[2]
                    if suffix:
                        fancy = '{0!s}_{1!s}'.format(fancy, suffix)
                    tools[fancy] = tool_factory(fancy, name, driver)
        except (subprocess.CalledProcessError, OSError):
            pass

    if not tools:
        errmsg = "Failed to load v5 tools"
        logger.debug(errmsg)
        raise GromacsToolLoadingError(errmsg)
    logger.debug("Loaded {0} v5 tools successfully!".format(len(tools)))
    return tools


def load_v4_tools():
    """ Load Gromacs 4.x tools automatically using some heuristic.

    Tries to load tools (1) in configured tool groups (2) and fails back  to
    automatic detection from ``GMXBIN`` (3) then to a prefilled list.

    Also load any extra tool configured in ``~/.gromacswrapper.cfg``

    :return: dict mapping tool names to GromacsCommand classes
    """
    logger.debug("Loading v4 tools...")

    names = config.get_tool_names()

    if len(names) == 0 and 'GMXBIN' in os.environ:
        names = find_executables(os.environ['GMXBIN'])

    if len(names) == 0 or len(names) > len(V4TOOLS) * 4:
        names = list(V4TOOLS)

    names.extend(config.get_extra_tool_names())

    tools = {}
    for name in names:
        fancy = make_valid_identifier(name)
        tools[fancy] = tool_factory(fancy, name, None)

    if not tools:
        errmsg = "Failed to load v4 tools"
        logger.debug(errmsg)
        raise GromacsToolLoadingError(errmsg)
    logger.debug("Loaded {0} v4 tools successfully!".format(len(tools)))
    return tools


def merge_ndx(*args):
    """ Takes one or more index files and optionally one structure file and
    returns a path for a new merged index file.

    :param args: index files and zero or one structure file
    :return: path for the new merged index file
    """
    ndxs = []
    struct = None
    for fname in args:
        if fname.endswith('.ndx'):
            ndxs.append(fname)
        else:
            if struct is not None:
                raise ValueError("only one structure file supported")
            struct = fname

    fd, multi_ndx = tempfile.mkstemp(suffix='.ndx', prefix='multi_')
    os.close(fd)
    atexit.register(os.unlink, multi_ndx)

    if struct:
        make_ndx = registry['Make_ndx'](f=struct, n=ndxs, o=multi_ndx)
    else:
        make_ndx = registry['Make_ndx'](n=ndxs, o=multi_ndx)

    _, _, _ = make_ndx(input=['q'], stdout=False, stderr=False)
    return multi_ndx


# Load tools
if config.MAJOR_RELEASE == '5':
    logger.debug("Trying to load configured Gromacs major release {0}".format(
        config.MAJOR_RELEASE))
    registry = load_v5_tools()
elif config.MAJOR_RELEASE == '4':
    logger.debug("Trying to load configured Gromacs major release {0}".format(
        config.MAJOR_RELEASE))
    registry = load_v4_tools()
else:
    logger.debug("No major release configured: trying 5 -> 4")
    try:
        registry = load_v5_tools()
    except GromacsToolLoadingError:
        try:
            registry = load_v4_tools()
        except GromacsToolLoadingError:
            errmsg = "Autoloading was unable to load any Gromacs tool"
            logger.critical(errmsg)
            raise GromacsToolLoadingError(errmsg)

# Aliases command names to run unmodified GromacsWrapper scripts on a machine
# with only 5.x
for fancy, cmd in registry.items():
    for c5, c4 in NAMES5TO4.iteritems():
        # have to check each one, since it's possible there are suffixes
        # like for double precision; cmd.command_name is Gromacs name
        # (e.g. 'convert-tpr') so we need to be careful in the processing below.
        name = cmd.command_name
        if name.startswith(c5):
            if c4 == c5:
                break
            else:
                # maintain suffix (note: need to split with fancy because Gromacs
                # names (c5) may contain '-' etc)
                name = c4 + fancy.split(make_valid_identifier(c5))[1]
                registry[make_valid_identifier(name)] = registry[fancy]
                break
    else:
        # the common case of just adding the 'g_'
        registry['G_{0!s}'.format(fancy.lower())] = registry[fancy]


# Patching up commands that may be useful to accept multiple index files
for name4, name5 in [('G_mindist', 'Mindist'), ('G_dist', 'Distance')]:
    if name4 in registry:
        cmd = registry[name4]
        registry[name4] = tool_factory(name4, cmd.command_name, cmd.driver,
                                       GromacsCommandMultiIndex)
        if name5 in registry:
            registry[name5] = registry[name4]


# Append class doc for each command
for name in registry.iterkeys():
    __doc__ += ".. class:: {0!s}\n    :noindex:\n".format(name)


# Finally add command classes to module's scope
globals().update(registry)
__all__ = ['GromacsCommandMultiIndex', 'merge_ndx']
__all__.extend(registry.keys())

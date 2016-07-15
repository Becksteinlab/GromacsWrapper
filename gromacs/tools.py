# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
:mod:`gromacs.tools` -- Gromacs commands classes
================================================

A Gromacs command class acts as a factory function that produces an
instance of a gromacs command (:class:`gromacs.core.GromacsCommand`)
with initial default values.

By convention, a class has the capitalized name of the corresponding Gromacs
tool; dots and dashes are replaced by underscores to make it a valid python
identifier. Gromacs 5 tools (e.g, `sasa`) are aliased to their Gromacs 4 tool
names (e.g, `g_sas`) for backwards compatibility.

The list of Gromacs tools to be loaded is configured with the ``tools`` and
``groups`` options of the ``~/.gromacswrapper.cfg`` file. If these options
are not provided guesses are made.

It is also possible to extend the basic commands and patch in additional
functionality. For example, the :class:`GromacsCommandMultiIndex` class makes a
command accept multiple index files and concatenates them on the fly; the
behaviour mimics Gromacs' "multi-file" input that has not yet been enabled for
all tools.

.. autoclass:: GromacsCommandMultiIndex

Example
-------

In this example we create two instances of the :class:`gromacs.tools.Trjconv`
command (which runs the Gromacs ``trjconv`` command)::

  import gromacs.tools as tools

  trjconv = tools.Trjconv()
  trjconv_compact = tools.Trjconv(ur='compact', center=True, boxcenter='tric', pbc='mol',
                                  input=('protein','system'),
                                  doc="Returns a compact representation of the system centered on the protein")

The first one, ``trjconv``, behaves as the standard commandline tool but the
second one, ``trjconv_compact``, will by default create a compact
representation of the input data by taking into account the shape of the unit
cell. Of course, the same effect can be obtained by providing the corresponding
arguments to ``trjconv`` but by naming the more specific command differently
one can easily build up a library of small tools that will solve a specifi,
repeatedly encountered problem reliably. This is particularly helpful when doing
interactive work.

Gromacs tools
-------------
.. The docs for the tool classes are auto generated.

.. autoclass:: Mdrun
   :members:
"""
from __future__ import absolute_import
__docformat__ = "restructuredtext en"

import os.path
import tempfile
import subprocess
import atexit

from . import config, exceptions
from .core import GromacsCommand


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


ALIASES5TO4 = {
    'grompp': 'grompp',
    'eneconv': 'eneconv',
    'sasa': 'g_sas',
    'distance': 'g_dist',
    'convert_tpr': 'tpbconv',
    'editconf': 'editconf',
    'pdb2gmx': 'pdb2gmx',
    'trjcat': 'trjcat',
    'trjconv': 'trjconv',
    'trjorder': 'trjorder',
    'xpm2ps': 'xpm2ps',
    'mdrun': 'mdrun',
    'make_ndx': 'make_ndx',
    'make_edi': 'make_edi',
    'dump': 'gmxdump',
    'check': 'gmxcheck',
    'genrestr': 'genrestr',
    'genion': 'genion',
    'genconf': 'genconf',
    'do_dssp': 'do_dssp',
    'solvate': 'genbox',
}


class GromacsCommandMultiIndex(GromacsCommand):
        def __init__(self, **kwargs):
            kwargs = self._fake_multi_ndx(**kwargs)
            super(GromacsCommandMultiIndex, self).__init__(**kwargs)

        def run(self,*args,**kwargs):
            kwargs = self._fake_multi_ndx(**kwargs)
            return super(GromacsCommandMultiIndex, self).run(*args, **kwargs)

        def _fake_multi_ndx(self, **kwargs):
            ndx = kwargs.get('n')
            if not (ndx is None or type(ndx) is basestring):
                if len(ndx) > 1:
                    if 's' in kwargs:
                        ndx.append(kwargs.get('s'))
                    kwargs['n'] = merge_ndx(*ndx)
            return kwargs


def tool_factory(clsname, name, driver, doc=None):
    """ GromacsCommand derived type factory. """
    clsdict = {
        'command_name': name,
        'driver': driver
    }
    if doc:
        clsdict['__doc__'] = doc
    return type(clsname, (GromacsCommand,), clsdict)


def find_executables(path):
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
    fails back to automatic detection from ``GMXBIN`` (3) then to rough guesses.

    In all cases the command ``gmx help`` is ran to get all tools available.

    :return: dict mapping tool names to GromacsCommand classes
    """
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
                    fancy = name.replace('-', '_').capitalize()
                    suffix = driver.partition('_')[2]
                    if suffix:
                        fancy = '%s_%s' % (fancy, suffix)
                    tools[fancy] = tool_factory(fancy, name, driver)
        except (subprocess.CalledProcessError, OSError):
            pass

    if len(tools) == 0:
        raise exceptions.GromacsToolLoadingError("Failed to load v5 tools")
    return tools


def load_v4_tools():
    """ Load Gromacs 4.x tools automatically using some heuristic.

    Tries to load tools (1) in configured tool groups (2) and fails back  to
    automatic detection from ``GMXBIN`` (3) then to a prefilled list.

    :return: dict mapping tool names to GromacsCommand classes
    """
    names = config.get_tool_names()

    if len(names) == 0 and 'GMXBIN' in os.environ:
        names = find_executables(os.environ['GMXBIN'])

    if len(names) == 0 or len(names) > len(V4TOOLS) * 4:
        names = V4TOOLS[:]

    tools = {}
    for name in names:
        try:
            null = open(os.devnull, 'w')
            subprocess.check_call([name], stdout=null, stderr=null)
        except subprocess.CalledProcessError:
            pass
        fancy = name.capitalize().replace('.', '_')
        tools[fancy] = tool_factory(fancy, name, None)

    if len(tools) == 0:
        raise exceptions.GromacsToolLoadingError("Failed to load v4 tools")
    return tools


def load_extra_tools():
    """ Load extra tools.

    :return: dict mapping tool names to GromacsCommand classes
    """
    names = config.cfg.get("Gromacs", "extra").split()
    tools = {}
    for name in names:
        fancy = name.capitalize().replace('-', '_')
        driver = 'gmx' if config.MAJOR_RELEASE == '5' else None
        tools[name] = tool_factory(fancy, name, driver)
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
            assert struct is None, "only one structure file supported"
            struct = fname

    fd, multi_ndx = tempfile.mkstemp(suffix='.ndx', prefix='multi_')
    atexit.register(os.unlink, multi_ndx)

    if struct:
        make_ndx = registry['Make_ndx'](f=struct, n=ndxs, o=multi_ndx)
    else:
        make_ndx = registry['Make_ndx'](n=ndxs, o=multi_ndx)

    _, _, _ = make_ndx(input=['q'], stdout=False, stderr=False)
    return multi_ndx


if config.MAJOR_RELEASE == '5':
    registry = load_v5_tools()
elif config.MAJOR_RELEASE == '4':
    registry = load_v4_tools()
else:
    try:
        registry = load_v5_tools()
    except exceptions.GromacsToolLoadingError:
        try:
            registry = load_v4_tools()
        except exceptions.GromacsToolLoadingError:
            raise exceptions.GromacsToolLoadingError("Unable to load any tool")

registry.update(load_extra_tools())


# Aliases command names to run unmodified GromacsWrapper scripts on a machine
# without Gromacs 4.x
for name in registry.copy():
    for c4, c5 in ALIASES5TO4.iteritems():
        # have to check each one, since it's possible there are suffixes
        # like for double precision
        if name.startswith(c5):
            # mantain suffix
            registry[c4 + name.split(c5)] = registry[name]
            break
    else:
        # the common case of just adding the 'g_'
        registry['G_%s' % name.lower()] = registry[name]

for name4, name5 in [('G_mindist', 'Mindist'), ('G_dist', 'Distance')]:
    if name4 in registry:
        klass = type(name4, (GromacsCommandMultiIndex,),  {
            'command_name': registry[name5].command_name,
            'driver': registry[name5].driver,
            '__doc__': registry[name5].__doc__
        })
        registry[name4] = klass
        if config.MAJOR_RELEASE == '5':
            registry[name5] = klass


# 5.0.5 compatibility hack
if 'Convert_tpr' in registry:
    registry['Tpbconv'] = registry['Convert_tpr']

# Append class doc for each command
for cmd in registry.itervalues():
    __doc__ += ".. class:: %s\n    :noindex:\n" % cmd.__name__


# finally add command classes to module's scope
globals().update(registry)
__all__ = ['GromacsCommandMultiIndex', 'merge_ndx']
__all__ = __all__ + registry.keys()

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
tool; dots are replaced by underscores to make it a valid python identifier.
Gromacs 5 tools (e.g, `sasa`) are aliased to their Gromacs 4 tool names
(e.g, `g_sas`) for backwards compatibility.

The list of Gromacs tools to be loaded is configured in
:data:`gromacs.config.gmx_tool_groups`.

Example
-------

In this example we create two instances of the :class:`gromacs.tools.Trjconv`
command (which runs the Gromacs ``trjconv`` command)::

  import gromacs.tools as tools

  trjconv = tools.Trjconv()
  trjconv_compact = tools.Trjconv(ur='compact', center=True, boxcenter='tric', pbc='mol',
                                  input=('protein','system'),
                                  doc="Returns a compact representation of the"
                                      "system centered on the protein")

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

import os.path, tempfile, subprocess, atexit, warnings

from . import config, exceptions
from .core import GromacsCommand


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


TOOLS_V4 = ("do_dssp", "editconf", "eneconv", "g_anadock", "g_anaeig",
            "g_analyze", "g_angle", "g_bar", "g_bond", "g_bundle", "g_chi",
            "g_cluster", "g_clustsize", "g_confrms", "g_covar", "g_current",
            "g_density", "g_densmap", "g_densorder", "g_dielectric",
            "g_dipoles", "g_disre", "g_dist", "g_dos", "g_dyecoupl", "g_dyndom",
            "genbox", "genconf", "g_enemat", "g_energy", "genion", "genrestr",
            "g_filter", "g_gyrate", "g_h2order", "g_hbond", "g_helix",
            "g_helixorient", "g_hydorder", "g_kinetics", "g_lie", "g_luck",
            "g_mdmat", "g_membed", "g_mindist", "g_morph", "g_msd",
            "gmxcheck", "gmxdump", "g_nmeig", "g_nmens", "g_nmtraj", "g_options",
            "g_order", "g_pme_error", "g_polystat", "g_potential",
            "g_principal", "g_protonate", "g_rama", "g_rdf", "g_rms",
            "g_rmsdist", "g_rmsf", "grompp", "g_rotacf", "g_rotmat",
            "g_saltbr", "g_sans", "g_sas", "g_select", "g_sgangle", "g_sham",
            "g_sigeps", "g_sorient", "g_spatial", "g_spol", "g_tcaf",
            "g_traj", "g_tune_pme", "g_vanhove", "g_velacc", "g_wham",
            "g_wheel", "g_x2top", "g_xrama", "make_edi", "make_ndx", "mdrun",
            "mk_angndx", "ngmx", "pdb2gmx", "tpbconv", "trjcat", "trjconv",
            "trjorder", "xpm2p")


version = config.cfg.get('Gromacs', 'release')
major_release = version.split('.')[0]


def tool_factory(clsname, name, driver):
    """ GromacsCommand derived class factory. """
    return type(clsname, (GromacsCommand,), {
        'command_name':name,
        'driver': driver
    })


def append_suffix(name):
    """ Append (or not) the optional command suffix. """
    suffix = config.cfg.get('Gromacs', 'suffix')
    if suffix:
        name += '_' + suffix
    return name


def load_v5_tools():
    """ Load Gromacs 5.x tools.
    :return: dict mapping tool names to GromacsCommand classes
    """
    driver = append_suffix('gmx')
    try:
        out = subprocess.check_output([driver, '-quiet', 'help', 'commands'])
    except subprocess.CalledProcessError:
        raise exceptions.GromacsToolLoadingError("Failed to load v5 tools")

    tools = {}
    for line in str(out).encode('ascii').splitlines()[5:-1]:
        if line[4] != ' ':
            name = line[4:line.index(' ', 4)]
            fancy = name.replace('-', '_').capitalize()
            tools[fancy] = tool_factory(fancy, name, driver)
    return tools


def load_v4_tools():
    """ Load Gromacs 4.x tools.
    :return: dict mapping tool names to GromacsCommand classes
    """
    names = config.cfg.get("Gromacs", "tools")
    if not names:
        names = [append_suffix(t) for t in TOOLS_V4]
    else:
        names = names.split()

    tools = {}
    for tool in names:
        fancy = tool.capitalize()
        tools[fancy] = tool_factory(fancy, tool, None)

    try:
        null = open(os.devnull, 'w')
        subprocess.check_call(['g_luck'], stdout=null, stderr=null)
    except subprocess.CalledProcessError:
        raise exceptions.GromacsToolLoadingError("Failed to load v4 tools")
    return tools


def load_extra_tools():
    """ Load extra tools.
    :return: dict mapping tool names to GromacsCommand classes
    """
    names = config.cfg.get("Gromacs", "extra").split()
    driver = append_suffix('gmx') if major_release == '5' else None
    tools = {}
    for name in names:
        fancy = name.capitalize().replace('-', '_')
        tools[name] = tool_factory(fancy, name, driver)
    return tools


if major_release == '5':
    registry = load_v5_tools()
elif major_release == '4':
    registry = load_v4_tools()
else:
    raise exceptions.GromacsToolLoadingError("Unknow Gromacs version %s" %
                                           version)
registry.update(load_extra_tools())


# Append class doc for each command
for cmd in registry.itervalues():
    __doc__ += ".. class:: %s\n    :noindex:\n" % cmd.__name__


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


class GromacsCommandMultiIndex(GromacsCommand):
        def __init__(self, **kwargs):
            kwargs = self._fake_multi_ndx(**kwargs)
            super(GromacsCommandMultiIndex, self).__init__(**kwargs)

        def run(self,*args,**kwargs):
            kwargs = self._fake_multi_ndx(**kwargs)
            return super(GromacsCommandMultiIndex, self).run(*args, **kwargs)

        def _fake_multi_ndx(self, **kwargs):
            ndx = kwargs.get('n')
            if not (ndx is None or type(ndx) is str):
                if len(ndx) > 1:
                    kwargs['n'] = merge_ndx(ndx, kwargs.get('s'))
            return kwargs


for name4, name5 in [('G_mindist', 'Mindist'), ('G_dist', 'Distance')]:
    if name5 in registry:
        klass = type(name4, (GromacsCommandMultiIndex,),  {
            'command_name': registry[name5].command_name,
            'driver': registry[name5].driver,
            '__doc__': registry[name5].__doc__
        })
        registry[name4] = klass
        registry[name5] = klass


# 5.0.5 compatibility hack
if 'Convert_tpr' in registry:
    registry['Tpbconv'] = registry['Convert_tpr']


# finally add command classes to module's scope
globals().update(registry)
__all__ = registry.keys()
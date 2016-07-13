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
Gromacs 5 tools (e.g, `sasa`) are aliased to their Gromacs 4 tool names (e.g, `g_sas`)
for backwards compatibility.

The list of Gromacs tools to be loaded is configured in
:data:`gromacs.config.gmx_tool_groups`.

It is also possible to extend the basic commands and patch in additional
functionality. For example, the :class:`GromacsCommandMultiIndex` class makes a
command accept multiple index files and concatenates them on the fly; the
behaviour mimics Gromacs' "multi-file" input that has not yet been enabled for
all tools.

.. autoclass:: GromacsCommandMultiIndex
   :members: run, _fake_multi_ndx, __del__

Example
-------

In this example we create two instances of the :class:`gromacs.tools.Trjconv` command (which
runs the Gromacs ``trjconv`` command)::

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
    'gmxdump': 'gmxdump',
    'gmxcheck': 'gmxcheck',
    'genrestr': 'genrestr',
    'genion': 'genion',
    'genconf': 'genconf',
    'do_dssp': 'do_dssp',
    'solvate': 'genbox',
}


def tool_factory(clsname, name, driver, doc=None):
    """ GromacsCommand derived type factory. """
    clsdict = {
        'command_name': name,
        'driver': driver
    }
    if doc:
        clsdict['__doc__'] = doc
    return type(clsname, (GromacsCommand,), clsdict)


def load_v5_tools():
    """ Load Gromacs 5.x tools automatically inferred from running ``gmx help``.

    :return: dict mapping tool names to GromacsCommand classes
    """
    try:
        out = subprocess.check_output(['gmx', '-quiet', 'help', 'commands'])
    except subprocess.CalledProcessError:
        raise exceptions.GromacsToolLoadingError("Failed to load v5 tools")

    tools = {}
    for line in str(out).encode('ascii').splitlines()[5:-1]:
        if line[4] != ' ':
            name = line[4:line.index(' ', 4)]
            fancy = name.replace('-', '_').capitalize()
            tools[fancy] = tool_factory(fancy, name, 'gmx')
    return tools


def load_v4_tools():
    """ Load Gromacs 4.x tools.

    :return: dict mapping tool names to GromacsCommand classes
    """
    tools = {}
    for name in config.get_tools():
        fancy = name.capitalize().replace('-', '_')
        tools[fancy] = tool_factory(fancy, name, None)
    try:
        null = open(os.devnull, 'w')
        subprocess.check_call(['grompp'], stdout=null, stderr=null)
    except subprocess.CalledProcessError:
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


if config.MAJOR_RELEASE == '5':
    registry = load_v5_tools()
elif config.MAJOR_RELEASE == '4':
    registry = load_v4_tools()
else:
    raise exceptions.GromacsToolLoadingError("Unknow Gromacs version %s" %
                                           config.RELEASE)
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

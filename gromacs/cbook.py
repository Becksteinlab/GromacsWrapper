# $Id$
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)

"""\
Gromacs Cook Book (cbook)
=========================

The ``cbook `` module contains short recipes for tasks that are often
repeated. In the simplest case this is just one of the gromacs tools
with a certain set of command line options.

By collecting these invocations here it, errors can be avoided by
repeatedly re-inventing the wheel and it can also be used as a
tutorial.
"""

# Right now the simplest thing to do is to just create instances with pre-set
# values; this works fine and is succinct but has some disadvantages:
# * the underlying gromacs tool is executed to extract the help string; this
#   adds to the import time
# * adding documentation is awkward
# 
# For more complicated cases one is probably better off by properly deriving a
# new class and set arguments explicitly in init (using kwargs['flag'] =
# default) ... or I can write some meta(??) class to do this nicely

import re

import gromacs
import tools

trj_compact = tools.G_trjconv(ur='compact', center=True, boxcenter='tric', pbc='mol',
                              input=('protein','system'),
                              doc="Returns a compact representation of the system centered on the protein")

def grompp_qtot(*args, **kwargs):
    """Run ``gromacs.grompp`` and return the total charge of the system."""

    # match '  System has non-zero total charge: -4.000001e+00',
    qtot_pattern = re.compile('System has non-zero total charge: *(?P<qtot>[-+]?\d*\.\d+[eE][-+]\d+)')
    kwargs['stdout'] = False
    rc, output, junk = gromacs.grompp(*args, **kwargs)
    qtot = 0
    for line in output.split('\n'):
        m = qtot_pattern.search(line)
        if m:
            qtot = float(m.group('qtot'))
            break
    return qtot


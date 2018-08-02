# GromacsWrapper: test_example.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

from __future__ import division, absolute_import, print_function

import contextlib
import os
import errno

import pytest

import gromacs.config
import gromacs.utilities

def get_GMXRC():
    # get GMXRC from installed Gromacs conda package
    for gmxexe in ('gmx', 'gmx_d', 'gmx_mpi', 'gmx_mpi_d', 'grompp', 'mdrun'):
        path = gromacs.utilities.which(gmxexe)
        if path is not None:
            break
    else:
        raise RuntimeError("Cannot find Gromacs installation")

    bindir = os.path.dirname(path)
    GMXRC = os.path.join(bindir, 'GMXRC')
    if not os.path.exists(GMXRC):
        raise IOError(errno.EEXIST, "Could not find Gromacs setup file", GMXRC)
    return GMXRC


@contextlib.contextmanager
def temp_environ():
    _environ = os.environ.copy()
    try:
        yield os.environ
    finally:
        os.environ.clear()
        os.environ.update(_environ)


def test_set_gmxrc_environment():
    # not threadsafe: function modifies the global process environment

    gmx_envvars = ('GMXBIN', 'GMXLDLIB', 'GMXMAN', 'GMXDATA',
                   'GMXPREFIX', 'GROMACS_DIR')

    GMXRC = get_GMXRC()

    with temp_environ() as environ:
        # clean environment so that we can detect changes
        for envvar in gmx_envvars:
            try:
                del environ[envvar]
            except KeyError:
                pass

        before = environ.copy()
        gromacs.config.set_gmxrc_environment(GMXRC)
        newvars = set(environ) - set(before)
        for envvar in gmx_envvars:
            assert envvar in newvars, \
                "GMX environment variable was not added correctly"

def test_check_setup():
    rc = gromacs.config.check_setup()
    assert rc in (True, False)

def test_get_configuration():
    cfg = gromacs.config.get_configuration()
    # could test more variables
    assert cfg.getpath('DEFAULT', 'configdir')


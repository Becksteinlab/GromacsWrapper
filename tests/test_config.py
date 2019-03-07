# GromacsWrapper: test_example.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

from __future__ import division, absolute_import, print_function

import contextlib
import os
import errno
import sys

if sys.version_info[0] < 3:
    from ConfigParser import NoOptionError, NoSectionError
else:
    from configparser import NoOptionError, NoSectionError

import pytest

import gromacs.config
import gromacs.utilities


@pytest.fixture
def GMXRC():
    # Try using GMXRC in config file:
    GMXRC = gromacs.config.cfg.get('Gromacs', 'gmxrc')
    if GMXRC:
        return GMXRC
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
        raise IOError(errno.ENOENT, "Could not find Gromacs setup file", GMXRC)
    return GMXRC


@contextlib.contextmanager
def temp_environ():
    _environ = os.environ.copy()
    try:
        yield os.environ
    finally:
        os.environ.clear()
        os.environ.update(_environ)


def test_set_gmxrc_environment(GMXRC):
    # not threadsafe: function modifies the global process environment

    gmx_envvars = ('GMXBIN', 'GMXLDLIB', 'GMXMAN', 'GMXDATA',
                   'GMXPREFIX', 'GROMACS_DIR')

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
    assert isinstance(cfg.getboolean('Gromacs', 'append_suffix'), bool)


def test_modified_config(modified_config):
    tools, append_suffix, Path = modified_config
    if tools != '':
        assert Path('~/gmx_mpi').expanduser().exists()
    assert gromacs.config.cfg.get('Gromacs', 'tools') == tools
    assert gromacs.config.cfg.get('Gromacs', 'append_suffix') == append_suffix


def test_get_boolean():
    # The code for getboolean was custom implemented based on the Python 3.7
    # ConfigParser.getboolean code.
    # These tests should be unnecessary for the Python 3 version of the code.
    cfg = gromacs.config.cfg
    assert isinstance(cfg.getboolean('Gromacs', 'append_suffix'), bool)
    assert isinstance(cfg.getboolean('Gromacs', 'append_suffix',
                                     fallback=True), bool)
    with pytest.raises(ValueError):
        cfg.getboolean('DEFAULT', 'configdir')
    with pytest.raises(ValueError):
        cfg.getboolean('DEFAULT', 'configdir', fallback=True)
    with pytest.raises(NoOptionError):
        cfg.getboolean('Gromacs', 'bool')
    with pytest.raises(NoSectionError):
        cfg.getboolean('Not a section', 'bool')
    cfg.set('Gromacs', 'bool', '')
    cfg.remove_option('Gromacs', 'bool')
    assert cfg.getboolean('Gromacs', 'bool', fallback=True) is True
    with pytest.raises(NoOptionError):
        cfg.getboolean('Gromacs', 'bool')

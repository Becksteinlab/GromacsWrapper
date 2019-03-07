# GromacsWrapper: test_example.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

# Ian Kenney's pytest-gmx plugin
# https://github.com/ianmkenney/pytest-gmx/

import distutils.spawn
import os
import pytest
import shutil
import sys

if sys.version_info[0] < 3:
    from ConfigParser import ConfigParser
    from pathlib2 import Path
else:
    from configparser import ConfigParser
    from pathlib import Path


def pytest_addoption(parser):
    """Add options to control Gromacs performance settings"""

    group = parser.getgroup('gmx')
    group.addoption(
            "--low-performance",
            action="store_true",
            dest='low_performance',
            help='Instruct Gromacs to run in low performance '
                 'mode (as defined in tests)',
    )
    group.addoption(
        '--no-append-suffix',
        action='store_false',
        dest='append_suffix',
        help='modify config file to set append_suffix to "no"'
    )
    group.addoption(
        '--link-gmx-mpi',
        action='store_true',
        dest='link_gmx_mpi',
        help='link the gmx executable to the home directory as gmx_mpi and '
             'add this path as a "tools" in the config file'
    )


@pytest.fixture(scope="session")
def low_performance(request):
    return request.config.option.low_performance


def gmx_mpi_linked(link):
    gmx_exe = distutils.spawn.find_executable('gmx')
    gmx_mpi = Path('~/gmx_mpi').expanduser()
    if not link:
        return ''
    else:
        gmx_mpi.symlink_to(gmx_exe)
        return str(gmx_mpi.expanduser())


@pytest.fixture
def modified_config(request):
    link_gmx_mpi = request.config.getoption('link_gmx_mpi')
    tools = str(Path('~/gmx_mpi').expanduser()) if link_gmx_mpi else ''
    append_suffix = 'yes' if request.config.getoption('append_suffix') else 'no'
    return tools, append_suffix, Path


path_config = Path('~/.gromacswrapper.cfg').expanduser()
gw_config = ConfigParser()
if path_config.exists():
    gw_config.read(str(path_config.resolve()))
    config_existed = True
else:
    gw_config.read('gromacs/templates/gromacswrapper.cfg')
    config_existed = False
config_backup = path_config.with_suffix('.bak')


def pytest_configure(config):
    link_gmx_mpi = config.getoption('link_gmx_mpi')
    append_suffix = 'yes' if config.getoption('append_suffix') else 'no'
    if config_existed:
        shutil.copy(str(path_config), str(config_backup))
    tools = gmx_mpi_linked(link_gmx_mpi)
    gw_config.set('Gromacs', 'tools', tools)
    gw_config.set('Gromacs', 'append_suffix', append_suffix)
    with open(str(path_config), 'w') as config_file:
        gw_config.write(config_file)


def pytest_unconfigure(config):
    if config_existed:
        config_backup.rename(str(path_config))
    else:
        os.remove(str(path_config))
    if config.option.link_gmx_mpi:
        gmx_mpi = Path('~/gmx_mpi').expanduser()
        gmx_mpi.unlink()

# GromacsWrapper: test_example.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

# Ian Kenney's pytest-gmx plugin
# https://github.com/ianmkenney/pytest-gmx/

import pytest

def pytest_addoption(parser):
    """Add options to control Gromacs performance settings"""

    group = parser.getgroup('gmx')
    group.addoption(
            "--low-performance",
            action="store_true",
            dest='low_performance',
            help='Instruct Gromacs to run in low performance mode (as defined in tests)',
    )

@pytest.fixture(scope="session")
def low_performance(request):
    return request.config.option.low_performance

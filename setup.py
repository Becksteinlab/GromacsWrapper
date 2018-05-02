# setuptools installation of GromacsWrapper
# Copyright (c) 2008-2011 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
#
# See the files INSTALL and README for details or visit
# https://github.com/Becksteinlab/GromacsWrapper
from __future__ import with_statement
from setuptools import setup, find_packages

import versioneer

with open("README.rst") as readme:
    long_description = readme.read()


setup(name="GromacsWrapper",
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      description="A python wrapper around the Gromacs tools.",
      long_description=long_description,
      author="Oliver Beckstein",
      author_email="orbeckst@gmail.com",
      license="GPLv3",
      url="https://github.com/Becksteinlab/GromacsWrapper",
      download_url="https://github.com/Becksteinlab/GromacsWrapper/downloads",
      keywords="science Gromacs analysis 'molecular dynamics'",
      classifiers=['Development Status :: 4 - Beta',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: GNU General Public License (GPL)',
                   'Operating System :: POSIX',
                   'Operating System :: MacOS :: MacOS X',
                   'Programming Language :: Python',
                   'Topic :: Scientific/Engineering :: Bio-Informatics',
                   'Topic :: Scientific/Engineering :: Chemistry',
                   'Topic :: Software Development :: Libraries :: Python Modules',
                   ],
      packages=find_packages(exclude=['tests','scripts','extras','doc/examples']),
      scripts = [
                 'scripts/gw-join_parts.py',
                 'scripts/gw-merge_topologies.py',
                 'scripts/gw-forcefield.py',
                 'scripts/gw-partial_tempering.py',
                 ],
      package_data={'gromacs': ['templates/*.sge', 'templates/*.pbs',  # template files
                                'templates/*.ll', 'templates/*.sh',
                                'templates/*.mdp', 'templates/*.cfg',
                                'tests/data/fileformats/top/*.mdp',    # test data
                                'tests/data/fileformats/top/*/*.top',
                                'tests/data/fileformats/top/*/*.gro',
                                'tests/data/*.log',
                                ],
                    },
      install_requires = ['numpy>=1.0',
                          'six',          # towards py 3 compatibility
                          'numkit',       # numerical helpers
                          ],
      # also currently requires
      #
      #    pip install git+git://github.com/ianmkenney/pytest-gmx.git
      #
      # for the low_performance pytest fixture: allows running
      #
      #    pytest --low-performance
      #
      # so that test are run with 'gmx mdrun -nt 2' as necessary on travis
      # (instead of the maximum number of threads, as determined by Gromacs)
      tests_require = ['pytest', 'numpy>=1.0', 'pandas>=0.17'],
      zip_safe = True,
)



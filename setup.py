# setuptools installation of GromacsWrapper
# Copyright (c) 2008-2011 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
#
# See the files INSTALL and README for details or visit
# https://github.com/Becksteinlab/GromacsWrapper
from __future__ import with_statement

from setuptools import setup, find_packages

with open("README.rst") as readme:
    long_description = readme.read()

# Dynamically calculate the version based on gromacs.VERSION.
# (but requires that we can actually import the package BEFORE it is
# properly installed!)
version = __import__('gromacs.version').get_version()

setup(name="GromacsWrapper",
      version=version,
      description="A python wrapper around the gromacs tools.",
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
      scripts = ['scripts/gw-fit_strip_trajectories.py',
                 'scripts/gw-join_parts.py',
                 'scripts/gw-merge_topologies.py',
                 'scripts/gw-forcefield.py',
                 'scripts/gw-partial_tempering.py',
                 ],
      package_data={'gromacs': ['templates/*.sge', 'templates/*.pbs',  # template files
                                'templates/*.ll', 'templates/*.sh',
                                'templates/*.mdp', 'templates/*.cfg',
                                'tests/test_fileformats/test_top/data/*/*.top',
                                ],
                    'vmd': ['*.tcl'],                                  # server start in VMD
                    },
      install_requires = ['numpy>=1.0',
                          'scipy',        # numkit needs it
                          ],              # basic package (w/o analysis)
      extras_require = {
                'analysis': ['matplotlib>=0.91.3',
                             'RecSQL>=0.7',
                             ],
                'numkit': ['scipy'],
                },
      zip_safe = True,
)



# setuptools installation of GromacsWrapper
# Copyright (c) 2008-2011 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
#
# See the files INSTALL and README for details or visit
# https://github.com/Becksteinlab/GromacsWrapper
from setuptools import setup, find_packages

import versioneer

with open("README.rst") as readme:
    long_description = readme.read()


setup(
    name="GromacsWrapper",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="A Python wrapper around the GROMACS tools.",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    author="Oliver Beckstein",
    author_email="orbeckst@gmail.com",
    license="GPLv3",
    url="https://github.com/Becksteinlab/GromacsWrapper",
    download_url="https://github.com/Becksteinlab/GromacsWrapper/downloads",
    keywords="science Gromacs analysis 'molecular dynamics'",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "License :: OSI Approved :: BSD License",
        "Operating System :: POSIX",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows ",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    packages=find_packages(
        exclude=["scripts", "tests", "tests.*", "extras", "doc/examples"]
    ),
    scripts=[
        "scripts/gw-join_parts.py",
        "scripts/gw-merge_topologies.py",
        "scripts/gw-forcefield.py",
        "scripts/gw-partial_tempering.py",
    ],
    package_data={
        "gromacs": [
            "templates/*.sge",
            "templates/*.pbs",  # template files
            "templates/*.ll",
            "templates/*.sh",
            "templates/*.mdp",
            "templates/*.cfg",
        ],
    },
    install_requires=[
        "matplotlib",
        "numkit",  # numerical helpers
        "numpy>=1.0",
        'setuptools; python_version >= "3.12"',  # contains 'pkg_resources' for versioneer
    ],
    tests_require=["pytest", "pandas>=0.17"],
    zip_safe=True,
)

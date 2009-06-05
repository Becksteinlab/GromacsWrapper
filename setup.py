# $Id$
# setuptools installation of GromacsWrapper
# Copyright (c) 2007-2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 2 (or higher, your choice)


from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup, find_packages

setup(name="GromacsWrapper",
      version="0.0.1",
      description="A python wrapper around the gromacs tools.",
      long_description="""\
A primitive wrapper around the Gromacs tools until we have proper
python bindings.
""",
      author="Oliver Beckstein",
      author_email="orbeckst@gmail.com",
      license="GPLv3",
      url="http://sbcb.bioch.ox.ac.uk/oliver/software", # not set up yet
      keywords="science Gromacs",
      packages=find_packages(exclude=['tests','extras','doc/examples']),
      install_requires=['numpy>=1.0',    # ... not now but sooner or later
                        ],
)

      

#!/bin/bash

SERVERDIR=/sansom/public_html/sbcb/oliver
PACKAGES=$SERVERDIR/download/Python
DOCS=$SERVERDIR/software/GromacsWrapper

distribution () {
  python setup.py sdist
  python setup.py bdist_egg
  rsync -v dist/* $PACKAGES
}

make_epydocs() {
  epydoc -v -o doc/epydoc --html --name=GromacsWrapper \
         --url=http://sbcb.bioch.ox.ac.uk/oliver/software/GromacsWrapper/ \
         gromacs gromacs/analysis/plugins/ vmd/
}

make_sphinx () {
  (cd doc/sphinx && make html)
}

docs () {
  make_epydocs 
  make_sphinx
  rsync -vrP doc/epydoc $DOCS
  rsync -vrP doc/sphinx/build/html $DOCS
}



distribution
docs

  

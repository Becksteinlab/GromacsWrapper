#!/bin/bash

PACKAGE=GromacsWrapper
EPYDOC_DIRS="gromacs gromacs/analysis/plugins/ vmd/ edPDB/"
PDF=GromacsWrapper.pdf

SERVERDIR=/sansom/public_html/html/sbcb/oliver
PACKAGES=$SERVERDIR/download/Python
DOCS=$SERVERDIR/software/$PACKAGE

COPY=1

usage="usage: $0 [OPTIONS] [cmd1 cmd2 ...]

Build distribution and python packages for $PACKAGE 
and copy files to the server directory (must be nfs mounted). 

default commands: 'distribution docs'

cmd           description
-----         -----------------
distribution  make sdist and egg, copy to $PACKAGES
dist
docs          'make_sphinx'
make_epydoc   source code docs, copy to $DOCS/epydoc # BROKEN
make_sphinx   documentation, copy to $DOCS/html


Options

-h           help
-c           publish on server (copy) [$COPY]
-n           do not copy              [$((1-COPY))]
-s DIR       server dir [${SERVERDIR}]
-p VERSION   python version, eg 2.5
"

function die () {
    echo 1>&2 "ERROR: $1"
    exit ${2:-1}
}

RSYNC () {
  if [ $COPY = 1 ]; then
      rsync $*;
  fi
}

distribution () {
  $PYTHON setup.py sdist \
      && $PYTHON setup.py bdist_egg \
      && RSYNC -v --checksum dist/* $PACKAGES \
      || die "Failed distribution"
}

dist () {
    distribution
}

make_epydocs() {
  epydoc -v -o doc/epydoc --html --name=$PACKAGE \
         --url=http://sbcb.bioch.ox.ac.uk/oliver/software/$PACKAGE/ \
         ${EPYDOC_DIRS}  \
      || die "Failed making epydoc"
  RSYNC -vrP --delete doc/epydoc $DOCS
}

make_sphinx () {
  (cd doc/sphinx && make clean && make html) || die "Failed making sphinx docs"
  echo "Created doc/html"
  RSYNC -vrP --delete doc/html $DOCS
}

make_pdf () {
  (cd doc/sphinx && make latex && cd build/latex && make all-pdf) || die "Failed making sphinx pdf"
  cp doc/sphinx/build/latex/$PDF doc
  echo "Updated pdf doc/$PDF"
}

sphinx () {
  make_sphinx
  make_pdf
}    

docs () {
  #make_epydocs 
  make_sphinx
  #make_pdf
}


COPY=1
while getopts hns:p: OPT; do
    case "$OPT" in
	h) echo "$usage"; exit 0;;
	n) COPY=0;;
	c) COPY=1;;
	s) SERVERDIR=$OPTARG;;
	p) PYVERSION=$OPTARG;;
	[?]) echo "Illegal option. See -h for usage.";
	     exit 1;;
    esac
done
shift $((OPTIND-1))


PACKAGES=$SERVERDIR/download/Python
DOCS=$SERVERDIR/software/$PACKAGE

case "$PYVERSION" in
   2.5|2.5.*)  PYTHON=python2.5;;
   2.6|2.6.*)  PYTHON=python2.6;;
   2.[0-4])    die "pyversion $PYVERSION not supported";;
   *)          PYTHON=python;;
esac

commands="$@"
[ -n "$commands" ] || commands="make_sphinx distribution"

for cmd in $commands; do
    eval "$cmd"
done


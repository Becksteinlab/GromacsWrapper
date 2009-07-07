#!/bin/bash

# XXX: would be nice to get version from setup.py and patch it into
# XXX: doc/sphinx/source/conf.py

PACKAGE=GromacsWrapper

SERVERDIR=/sansom/public_html/sbcb/oliver
PACKAGES=$SERVERDIR/download/Python
DOCS=$SERVERDIR/software/$PACKAGE

usage="usage: $0 [OPTIONS] [cmd1 cmd2 ...]

Build distribution and python packages for $PACKAGE 
and copy files to the server directory (must be nfs mounted). 

default commands: 'distribution docs'

cmd           description
-----         -----------------
distribution  make sdist and egg, copy to $PACKAGES
docs          'make_epydoc make_sphinx'
make_epydoc   source code docs, copy to $DOCS/epydoc
make_sphinx   documentation, copy to $DOCS/html


Options

-h           help
-n           do not copy
-s DIR       server dir [${SERVERDIR}]
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
  python setup.py sdist \
      && python setup.py bdist_egg \
      && RSYNC -v --checksum dist/* $PACKAGES \
      || die "Failed distribution"
}

make_epydocs() {
  epydoc -v -o doc/epydoc --html --name=$PACKAGE \
         --url=http://sbcb.bioch.ox.ac.uk/oliver/software/$PACKAGE/ \
         gromacs gromacs/analysis/plugins/ vmd/  \
      || die "Failed making epydoc"
  RSYNC -vrP --delete doc/epydoc $DOCS
}

make_sphinx () {
  (cd doc/sphinx && make html) || die "Failed making sphinx docs"
  RSYNC -vrP --delete doc/sphinx/build/html $DOCS
}

docs () {
  make_epydocs 
  make_sphinx
}


COPY=1
while getopts hns: OPT; do
    case "$OPT" in
	h) echo "$usage"; exit 0;;
	n) COPY=0;;
	s) SERVERDIR=$OPTARG;;
	[?]) echo "Illegal option. See -h for usage.";
	     exit 1;;
    esac
done
shift $((OPTIND-1))

PACKAGES=$SERVERDIR/download/Python
DOCS=$SERVERDIR/software/$PACKAGE


commands="$@"
[ -n "$commands" ] || commands="distribution docs"

for cmd in $commands; do
    eval "$cmd"
done


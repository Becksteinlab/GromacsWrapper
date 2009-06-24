#!/bin/sh
python setup.py sdist
python setup.py bdist_egg
rsync -v dist/* /sansom/public_html/sbcb/oliver/download/Python/

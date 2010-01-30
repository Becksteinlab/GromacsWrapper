# utilities for edPDB
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
:mod:`edPDB.utilities` -- Helper functions and classes
========================================================

The module defines some convenience functions and classes that are
used in other modules


Functions
---------

Functions that improve list processing and which do *not* treat
strings as lists:

.. autofunction:: iterable
.. autofunction:: asiterable

"""


def iterable(obj):
    """Returns ``True`` if *obj* can be iterated over and is *not* a  string."""
    if type(obj) is str:
        return False    # avoid iterating over characters of a string

    if hasattr(obj, 'next'):
        return True    # any iterator will do 
    try: 
        len(obj)       # anything else that might work
    except TypeError: 
        return False
    return True

def asiterable(obj):
    """Returns obj so that it can be iterated over; a string is *not* treated as iterable"""
    if not iterable(obj):
        obj = [obj]
    return obj

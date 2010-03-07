# GromacsWrapper: collections.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under GPL3+

"""
:mod:`analysis.collections` -- Handling of groups of simulation instances
=========================================================================

This module contains classes and functions that combine multiple
:class:`gromacs.analysis.core.Simulation` objects. In this way the
same kind of analysis or plotting task can be carried out
simultaneously for all simulations in the collection.

.. autoclass:: Collection
"""

import os.path
import cPickle
from numpy import all, any

class Collection(list):
    """Multiple objects (organized as a list).

    Methods are applied to all objects in the Collection and returned
    as new Collection:

      >>> from gromacs.analysis.collections import Collection
      >>> animals = Collection(['ant', 'boar', 'ape', 'gnu'])
      >>> animals.startswith('a')
      Collection([True, False, True, False])

    Similarly, attributes are returned as a Collection.

    Using :meth:`Collection.save` one can save the whole collection to
    disk and restore it later with the :meth:`Collection.load` method

      >>> animals.save('zoo')
      >>> arc = Collection()
      >>> arc.load('zoo')
      >>> arc.load('zoo', append=True)
      >>> arc
      ['ant', 'boar', 'ape', 'gnu', 'ant', 'boar', 'ape', 'gnu']
    """
    # note: do not use with multiple inheritance -- why, I think it could work now...

    def save(self, filename):
        """Pickle the whole collection to *filename*.
        
        If no extension is provided, ".collection" is appended.
        """
        cPickle.dump(self, open(self._canonicalize(filename), 'wb'), 
                     protocol=cPickle.HIGHEST_PROTOCOL)

    def load(self, filename, append=False):
        """Load collection from pickled file *filename*.

        *append* determines if the saved collection is added to the current one
        or if it replaces the current content.

        If no extension is provided, ".collection" is appended.
        """
        tmp = cPickle.load(open(self._canonicalize(filename), 'rb'))
        if append:
            self.extend(tmp)
        else:
            self[:] = tmp[:]
        del tmp

    def tolist(self):
        """Return contents as a simple list."""
        return self[:]

    def _canonicalize(self, filename):
        """Use .collection as extension unless provided"""
        path, ext = os.path.splitext(filename)
        if not ext:
            ext = ".collection"
        return path + ext

    def __getnewargs__(self, *args, **kwargs):
        """Provide proper initialization to make pickling with protocol 2 work"""
        return (self.tolist(),)

    def __getattribute__(self, attr):
        try:
            return super(Collection, self).__getattribute__(attr)
        except AttributeError:
            pass
        
        for o in self:
            failures = []
            if not hasattr(o, attr):
                failures.append(o)
        if len(failures) > 0:
            raise AttributeError("The following members of the collection do not "
                                 "implement the attribute %(attr)r:\n%(failures)r\n"
                                 % vars())

        # analyze attribute: functions (the ones with __call__) get delayed, simple
        # attributes are looked up immediately
        iscallable = [hasattr(o.__getattribute__(attr), '__call__') for o in self]
        if all(iscallable):
            def runall(*args, **kwargs):
                """Apply function to all members and return a new Collection"""
                return Collection([o.__getattribute__(attr)(*args, **kwargs) for o in self])
            runall.__name__ = attr
            return runall
        elif any(iscallable):
            raise TypeError("Attribute %r is callable only for some objects" % attr)

        return Collection([o.__getattribute__(attr) for o in self])

    def __add__(self, x):
        return Collection(super(Collection, self).__add__(x))

    def __repr__(self):
        return self.__class__.__name__+"(%s)" % super(Collection, self).__repr__()

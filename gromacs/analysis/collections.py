# $Id$
"""
:mod:`analysis.collections` -- Handling of groups of simulation instances
=========================================================================

This module contains classes and functions that combine multiple
:class:`gromacs.analysis.core.Simulation` objects. In this way the
same kind of analysis or plotting task can be carried out
simultaneously for all simulations in the collection.

.. autoclass:: Collection
"""

class Collection(list):
    """Multiple Simulation objects (organized as a list).

    Methods are applied to all simulation objects in the Collection.

    Attribute lookup does not work.
    """
    # note: do not use multiple inheritance because __init__ is iffy

    def __init__(self, *args):
        """Build the Collection from Simulation objects."""
        super(Collection, self).__init__(args) 

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
        def runall(*args, **kwargs):
            return [o.__getattribute__(attr)(*args, **kwargs) for o in self]                
        return runall
        

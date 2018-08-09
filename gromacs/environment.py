# GromacsWrapper environment.py
# Copyright (c) 2011 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
:mod:`gromacs.environment` -- Run time modification of behaviour
================================================================

Some aspects of GromacsWrapper can be determined globally. The
corresponding flags :class:`Flag` are set in the environment (think of
them like environment variables). They are accessible through the
pseudo-dictionary :data:`gromacs.environment.flags`.

The entries appear as 'name'-'value' pairs. Flags check values and illegal ones
raise a :exc:`ValueError`. Documentation on all flags can be obtained with ::

 print gromacs.environment.flags.doc()


List of GromacsWrapper flags with default values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: flagsDocs


Classes
~~~~~~~

.. data:: flags

.. autoclass:: Flags
   :members:

.. autoclass:: Flag
   :members:

"""
import six

# Flags infrastructure taken from MDAnalysis.core.__init__ (same author ... :-) )

# set up flags for core routines (more convoluted than strictly necessary but should
# be clean to add more flags if needed)
class Flags(dict):
    """Global registry of flags. Acts like a dict for item access.

    There are a number flags defined that influence how MDAnalysis behaves. They are
    accessible through the pseudo-dictionary

      :data:`gromacs.environment.flags`

    The entries appear as 'name'-'value' pairs. Flags check values and illegal ones
    raise a :exc:`ValueError`. Documentation on all flags can be obtained with ::

      print gromacs.environment.flags.__doc__

    New flags are added with the :meth:`Flags.register` method which takes a new :class:`Flag`
    instance as an argument.
    """
    def __init__(self,*args):
        """For **developers**: Initialize Flags registry with a *list* of :class:`Flag` instances."""
        super(Flags,self).__init__([(flag.name,flag) for flag in args])
    def get_flag(self,name):
        return super(Flags,self).__getitem__(name)
    def doc(self):
        """Shows doc strings for all flags."""
        return "\n\n".join([flag.__doc__ for flag in self._itervalues()])
    def register(self,flag):
        """Register a new :class:`Flag` instance with the Flags registry."""
        super(Flags,self).__setitem__(flag.name,flag)
    def update(self,*flags):
        """Update Flags registry with a list of :class:`Flag` instances."""
        super(Flags,self).update([(flag.name,flag) for flag in flags])
    def setdefault(self,k,d=None):
        raise NotImplementedError
    def __getitem__(self,name):
        return self.get_flag(name).get()
    def __setitem__(self,name,value):
        self.get_flag(name).set(value)
    def _itervalues(self):
        return six.itervalues(super(Flags,self))
    def _items(self):
        return super(Flags,self).items()
    def itervalues(self):
        for flag in self._itervalues():
            yield flag.value
    def iteritems(self):
        for flag in self._itervalues():
            yield flag.name,flag.value
    def values(self):
        return [flag.value for flag in self._itervalues()]
    def items(self):
        return [(flag.name,flag.value) for flag in self._itervalues()]
    def __repr__(self):
        return str(self.items())

class FlagsDynamicDocs(Flags):
    # docs are generated on the fly for interactive use; but because
    # this does not work well with the sphinx documentation system
    # ("AttributeError: 'property' object has no attribute
    # 'expandtabs'") we split the class...
    @property
    def __doc__(self):
        # generate dynamic docs on all flags
        return self.doc()


class IdentityMapping(dict):
    def __getitem__(self,key):
        return key

class Flag(object):
    """A Flag, essentially a variable that knows its default and legal values."""
    def __init__(self,name,default,mapping=None,doc=None):
        """Create a new flag which will be registered with FLags.

          newflag = Flag(name,default,mapping,doc)

        :Arguments:
         *name*
            name of the flag, must be a legal python name
         *default*
            default value
         *mapping*
            dict that maps allowed input values to canonical values;
            if ``None`` then no argument checking will be performed and
            all values are directly set.
         *doc*
            doc string; may contain string interpolation mappings for::

                    %%(name)s        name of the flag
                    %%(default)r     default value
                    %%(value)r       current value
                    %%(mapping)r     mapping

            Doc strings are generated dynamically and reflect the current state.
        """
        self.name = name
        self.value = default
        self.default = default
        # {v1:v1,v2:v1,v3:v3, ...} mapping of allowed values to canonical ones
        self.mapping = mapping or IdentityMapping()
        self._doctemplate = "**%(name)s** = *%(value)r*\n" + (doc or "*undocumented flag*")
    def get(self):
        return self.value
    def set(self,value):
        if value is not None:
            try:
                self.value = self.mapping[value]
            except KeyError:
                raise ValueError("flag must be None or one of "+str(self.mapping.keys()))
        return self.get()
    def prop(self):
        """Use this for property(**flag.prop())"""
        return {'fget':self.get, 'fset':self.set, 'doc':self.__doc__}
    def __repr__(self):
        return """Flag('{name!s}',{value!r})""".format(**self.__dict__)


class _Flag(Flag):
    @property
    def __doc__(self):
        # generate dynamic docs with current values
        return  self._doctemplate % self.__dict__

_flags = [
    _Flag('capture_output',
          False,
          {True: True,
           False: False,
           'file': 'file',
           },
          """
            Select if Gromacs command output is *always* captured.

            >>> flags['%(name)s'] = %(value)r

            By default a :class:`~gromacs.core.GromacsCommand` will
            direct STDOUT and STDERR output from the command itself to
            the screen (through /dev/stdout and /dev/stderr). When
            running the command, this can be changed with the keywords
            *stdout* and *stderr* as described in :mod:`gromacs.core`
            and :class:`~gromacs.core.Command`.

            If this flag is set to ``True`` then by default STDOUT and
            STDERR are captured as if one had set ::

               stdout=False, stderr=False

            Explicitly setting *stdout* and/or *stderr* overrides the
            behaviour described above.

            If set to the special keyword ``"file"` then the command
            writes to the file whose name is given by
            ``flags['capture_output_filename']``. This file is
            *over-written* for each command. In this way one can
            investigate the output from the last command (presumably
            because it failed). STDOUT and STDERR are captured into
            this file by default. STDERR is printed first and then
            STDOUT, which does not necessarily reflect the order of
            output one would see on the screen.

            The default is %(default)r.
          """
          ),
    _Flag('capture_output_filename',
          'gromacs_captured_output.txt',
          doc="""
            Name of the file that captures output if ``flags['capture_output'] = "file"

            >>> flags['%(name)s'] = %(value)r

            This is an *experimental* feature. The default is %(default)r.
          """),
    ]

#: Global flag registry for :mod:`gromacs.environment`.
#: Can be accessed like a dictionary and appears to the casual user as such.
flags = FlagsDynamicDocs(*_flags)
del _flags

# only for sphinx docs
class flagsDocs(object):
    __doc__ = flags.doc()

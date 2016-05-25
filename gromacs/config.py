# GromacsWrapper config.py
# Copyright (c) 2009-2011 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

""":mod:`gromacs.config` -- Configuration for GromacsWrapper
==========================================================

The config module provides configurable options for the whole package;
It mostly serves to define which gromacs tools and other scripts are
exposed in the :mod:`gromacs` package and where template files are
located. The user can configure *GromacsWrapper* by

1. editing the global configuration file ``~/.gromacswrapper.cfg``

2. placing template files into directories under ``~/.gromacswrapper``
   (:data:`gromacs.config.configdir`) which can be processed instead
   of the files that come with *GromacsWrapper*

In order to **set up a basic configuration file and the directories**
a user should execute :func:`gromacs.config.setup` at least once. It
will prepare the user configurable area in their home directory and it
will generate a default global configuration file
``~/.gromacswrapper.cfg`` (the name is defined in
:data:`CONFIGNAME`)::

  import gromacs
  gromacs.config.setup()

If the configuration file is edited then one can force a rereading of
the new config file with :func:`gromacs.config.get_configuration`::

 gromacs.config.get_configuration()

However, this will not update the available command classes (e.g. when new
executables were added to a tool group). In this case one either has to
:func:`reload` a number of modules (:mod:`gromacs`, :mod:`gromacs.config`,
:mod:`gromacs.tools`) although it is by far easier simply to quit python and
freshly ``import gromacs``.

Almost all aspects of *GromacsWrapper* (paths, names, what is loaded)
can be changed from within the configuration file. The only exception
is the name of the configuration file itself: This is hard-coded as
``~/.gromacswrapper.cfg`` although it is possible to read other
configuration files with the *filename* argument to
:func:`~gromacs.config.get_configuration`.


Configuration management
------------------------

Important configuration variables are

.. autodata:: configdir
.. autodata:: path

When GromacsWrapper starts up it runs :func:`check_setup`. This
notifies the user if any config files or directories are missing and
suggests to run :func:`setup`.  The check if the default set up exists
can be suppressed by setting the environment variable
:envvar:`GROMACSWRAPPER_SUPPRESS_SETUP_CHECK` to 'true' ('yes' and '1'
also work).


Users
~~~~~

Users will likely only need to run :func:`gromacs.config.setup` once and
perhaps occasionally execute :func:`gromacs.config.get_configuration`. Mainly
the user is expected to configure *GromacsWrapper* by editing the configuration
file ``~/.gromacswrapper.cfg`` (which has ini-file syntax as described in
:mod:`ConfigParser`).

.. autofunction:: setup
.. autofunction:: get_configuration
.. autofunction:: check_setup

Developers
~~~~~~~~~~

Developers are able to access all configuration data through
:data:`gromacs.config.cfg`, which represents the merger of the package default
values and the user configuration file values.

.. autodata:: cfg
.. autoclass:: GMXConfigParser
   :members:

A subset of important data is also made available as top-level package
variables as described under `Location of template files`_ (for historical
reasons); the same variable are also available in the dict
:data:`gromacs.config.configuration`.

.. autodata:: configuration

Default values are hard-coded in

.. autodata:: CONFIGNAME
.. autodata:: defaults

Accessing configuration and template files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following functions can be used to access configuration data. Note that
files are searched first with their full filename, then in all directories
listed in :data:`gromacs.config.path`, and finally within the package itself.

.. autofunction:: get_template
.. autofunction:: get_templates


Logging
-------

Gromacs commands log their invocation to a log file; typically at
loglevel *INFO* (see the python `logging module`_ for details).

.. autodata:: logfilename
.. autodata:: loglevel_console
.. autodata:: loglevel_file

.. _logging module: http://docs.python.org/library/logging.html


Gromacs tools and scripts
-------------------------

``load_*`` variables are lists that contain instructions to other
parts of the code which packages and scripts should be wrapped.

.. autodata:: load_tools

:data:`load_tools` is populated from lists of executable names in the
configuration file.

The tool groups are strings that contain white-space separated file
names of Gromacs tools. These lists determine which tools are made
available as classes in :mod:`gromacs.tools`.

A typical Gromacs tools section of the config file looks like this::

   [Gromacs]
   # Release of the Gromacs package to which information in this sections applies.
   release = 4.5.3

   # tools contains the file names of all Gromacs tools for which classes are generated.
   # Editing this list has only an effect when the package is reloaded.
   # (Note that this example has a much shorter list than the actual default.)
   tools =
         editconf make_ndx grompp genion genbox
         grompp pdb2gmx mdrun

   # Additional gromacs tools that should be made available.
   extra =
         g_count     g_flux
         a_gridcalc  a_ri3Dc  g_ri3Dc

   # which tool groups to make available as gromacs.NAME
   groups = tools extra

For Gromacs 5.x use a section like the following, where the driver
command ``gmx`` is added as a prefix::

   [Gromacs]
   # Release of the Gromacs package to which information in this sections applies.
   release = 5.0.5

   # tools contains the command names of all Gromacs tools for which classes are generated.
   # Editing this list has only an effect when the package is reloaded.
   # (Note that this example has a much shorter list than the actual default.)
   tools =
         gmx:editconf gmx:make_ndx gmx:grompp gmx:genion gmx:solvate
         gmx:insert-molecule gmx:convert-tpr
         gmx:grompp gmx:pdb2gmx gmx:mdrun

   # which tool groups to make available as gromacs.NAME
   groups = tools

For example, on the commandline you would run ::

   gmx grompp -f md.mdp -c system.gro -p topol.top -o md.tpr

and within GromacsWrapper this would become ::

   gromacs.grompp(f="md.mdp", c="system.gro", p="topol.top", o="md.tpr")

(The driver command is stripped and only the "command name" is used to
identify the command. This makes it easier to migrate GromacsWrapper
scripts from Gromacs 4.x to 5.x.).

The driver command can be changed per-tool, allowing for mpi
and non-mpi versions to be used. For example::

tools = gmx_mpi:mdrun gmx:pdb2gmx

.. Note:: Because of `changes in the Gromacs tool in 5.x`_,
          GromacsWrapper scripts might break, even if the tool
          names are still the same.

.. _`changes in the Gromacs tool in 5.x`:
   http://www.gromacs.org/Documentation/How-tos/Tool_Changes_for_5.0



Location of template files
--------------------------

*Template variables* list files in the package that can be used as
templates such as run input files. Because the package can be a zipped
egg we actually have to unwrap these files at this stage but this is
completely transparent to the user.

.. autodata:: qscriptdir
.. autodata:: templatesdir
.. autodata:: managerdir
.. autodata:: templates
.. autodata:: qscript_template

"""
from __future__ import absolute_import, with_statement

import os, errno
from ConfigParser import SafeConfigParser

from pkg_resources import resource_filename, resource_listdir

from . import utilities

# Defaults
# --------
# hard-coded package defaults

#: name of the global configuration file.
CONFIGNAME = os.path.expanduser(os.path.join("~",".gromacswrapper.cfg"))

#: Holds the default values for important file and directory locations.
#:
#: :data:`configdir`
#:    Directory to store user templates and configurations.
#:    The default value is ``~/.gromacswrapper``.
#: :data:`qscriptdir`
#:    Directory to store user supplied queuing system scripts as
#:    used by :mod:`gromacs.qsub`.
#:    The default value is ``~/.gromacswrapper/qscripts``.
#: :data:`templatesdir`
#:    Directory to store user supplied template files such as mdp files.
#:    The default value is ``~/.gromacswrapper/templates``.
#: :data:`managerdir`
#:    Directory to store configuration files for different queuing system
#:    managers as used in :mod:`gromacs.manager`.
#:    The default value is ``~/.gromacswrapper/managers``.
defaults = {
     'configdir': os.path.expanduser(os.path.join("~",".gromacswrapper")),
     'logfilename': "gromacs.log",
     'loglevel_console': 'INFO',
     'loglevel_file': 'DEBUG',
}
defaults['qscriptdir'] = os.path.join(defaults['configdir'], 'qscripts')
defaults['templatesdir'] = os.path.join(defaults['configdir'], 'templates')
defaults['managerdir'] = os.path.join(defaults['configdir'], 'managers')


# Logging
# -------
import logging
logger = logging.getLogger("gromacs.config")

#: File name for the log file; all gromacs command and many utility functions (e.g. in
#: :mod:`gromacs.cbook` and :mod:`gromacs.setup`) append messages there. Warnings and
#: errors are also recorded here. The default is *gromacs.log*.
logfilename = defaults['logfilename']

#: The default loglevel that is still printed to the console.
loglevel_console = logging.getLevelName(defaults['loglevel_console'])

#: The default loglevel that is still written to the :data:`logfilename`.
loglevel_file = logging.getLevelName(defaults['loglevel_file'])


# User-accessible configuration
# -----------------------------

# (written in a clumsy manner because of legacy code and because of
# the way I currently generate the documentation)

#: Directory to store user templates and rc files.
#: The default value is ``~/.gromacswrapper``.
configdir = defaults['configdir']

#: Directory to store user supplied queuing system scripts.
#: The default value is ``~/.gromacswrapper/qscripts``.
qscriptdir = defaults['qscriptdir']

#: Directory to store user supplied template files such as mdp files.
#: The default value is ``~/.gromacswrapper/templates``.
templatesdir = defaults['templatesdir']

#: Directory to store configuration files for remote queuing systems
#: :class:`gromacs.qsub.Manager` instances.
#: The default value is ``~/.gromacswrapper/managers``.
managerdir = defaults['managerdir']

#: List of all configuration directories.
config_directories = [configdir, qscriptdir, templatesdir, managerdir]


#: Search path for user queuing scripts and templates. The internal package-supplied
#: templates are always searched last via :func:`gromacs.config.get_templates`.
#: Modify :data:`gromacs.config.path` directly in order to customize the template
#: and qscript searching. By default it has the value ``['.', qscriptdir,
#: templatesdir]``.
#: (Note that it is not a good idea to have template files and qscripts with the
#: same name as they are both searched on the same path.)
#: :data:`path` is updated whenever cfg is re-read with :func:`get_configuration`.
path = [os.path.curdir, qscriptdir, templatesdir]


# Location of template files
# --------------------------

def _generate_template_dict(dirname):
    """Generate a list of included files *and* extract them to a temp space.

    Templates have to be extracted from the egg because they are used
    by external code. All template filenames are stored in
    :data:`config.templates`.
    """
    return dict((resource_basename(fn), resource_filename(__name__, dirname+'/'+fn))
                for fn in resource_listdir(__name__, dirname)
                if not fn.endswith('~'))

def resource_basename(resource):
    """Last component of a resource (which always uses '/' as sep)."""
    if resource.endswith('/'):
         resource = resource[:-1]
    parts = resource.split('/')
    return parts[-1]

templates = _generate_template_dict('templates')
"""*GromacsWrapper* comes with a number of templates for run input files
and queuing system scripts. They are provided as a convenience and
examples but **WITHOUT ANY GUARANTEE FOR CORRECTNESS OR SUITABILITY FOR
ANY PURPOSE**.

All template filenames are stored in
:data:`gromacs.config.templates`. Templates have to be extracted from
the GromacsWrapper python egg file because they are used by external
code: find the actual file locations from this variable.

**Gromacs mdp templates**

   These are supplied as examples and there is **NO GUARANTEE THAT THEY
   PRODUCE SENSIBLE OUTPUT** --- check for yourself!  Note that only
   existing parameter names can be modified with
   :func:`gromacs.cbook.edit_mdp` at the moment; if in doubt add the
   parameter with its gromacs default value (or empty values) and
   modify later with :func:`~gromacs.cbook.edit_mdp`.

   The safest bet is to use one of the ``mdout.mdp`` files produced by
   :func:`gromacs.grompp` as a template as this mdp contains all
   parameters that are legal in the current version of Gromacs.

**Queuing system templates**

   The queing system scripts are highly specific and you will need to add your
   own into :data:`gromacs.config.qscriptdir`.
   See :mod:`gromacs.qsub` for the format and how these files are processed.
"""

#: The default template for SGE/PBS run scripts.
qscript_template = templates['local.sh']


# Functions to access configuration data
# --------------------------------------

def get_template(t):
    """Find template file *t* and return its real path.

    *t* can be a single string or a list of strings. A string
    should be one of

    1. a relative or absolute path,
    2. a file in one of the directories listed in :data:`gromacs.config.path`,
    3. a filename in the package template directory (defined in the template dictionary
       :data:`gromacs.config.templates`) or
    4. a key into :data:`~gromacs.config.templates`.

    The first match (in this order) is returned. If the argument is a
    single string then a single string is returned, otherwise a list
    of strings.

    :Arguments: *t* : template file or key (string or list of strings)
    :Returns:   os.path.realpath(*t*) (or a list thereof)
    :Raises:    :exc:`ValueError` if no file can be located.

    """
    templates = [_get_template(s) for s in utilities.asiterable(t)]
    if len(templates) == 1:
         return templates[0]
    return templates

def get_templates(t):
    """Find template file(s) *t* and return their real paths.

    *t* can be a single string or a list of strings. A string should
    be one of

    1. a relative or absolute path,
    2. a file in one of the directories listed in :data:`gromacs.config.path`,
    3. a filename in the package template directory (defined in the template dictionary
       :data:`gromacs.config.templates`) or
    4. a key into :data:`~gromacs.config.templates`.

    The first match (in this order) is returned for each input argument.

    :Arguments: *t* : template file or key (string or list of strings)
    :Returns:   list of os.path.realpath(*t*)
    :Raises:    :exc:`ValueError` if no file can be located.

    """
    return [_get_template(s) for s in utilities.asiterable(t)]

def _get_template(t):
    """Return a single template *t*."""
    if os.path.exists(t):           # 1) Is it an accessible file?
         pass
    else:
         _t = t
         _t_found = False
         for d in path:              # 2) search config.path
              p = os.path.join(d, _t)
              if os.path.exists(p):
                   t = p
                   _t_found = True
                   break
         _t = os.path.basename(t)
         if not _t_found:            # 3) try template dirs
              for p in templates.values():
                   if _t == os.path.basename(p):
                        t = p
                        _t_found = True     # NOTE: in principle this could match multiple
                        break               #       times if more than one template dir existed.
         if not _t_found:            # 4) try it as a key into templates
              try:
                   t = templates[t]
              except KeyError:
                   pass
              else:
                   _t_found = True
         if not _t_found:            # 5) nothing else to try...
              raise ValueError("Failed to locate the template file %(t)r." % vars())
    return os.path.realpath(t)


class GMXConfigParser(SafeConfigParser):
     """Customized :class:`ConfigParser.SafeConfigParser`."""
     cfg_template = 'gromacswrapper.cfg'

     def __init__(self, *args, **kwargs):
          """Reads and parses the configuration file.

          Default values are loaded and then replaced with the values from
          ``~/.gromacswrapper.cfg`` if that file exists. The global
          configuration instance :data:`gromacswrapper.config.cfg` is updated
          as are a number of global variables such as :data:`configdir`,
          :data:`qscriptdir`, :data:`templatesdir`, :data:`logfilename`, ...

          Normally, the configuration is only loaded when the :mod:`gromacswrapper`
          package is imported but a re-reading of the configuration can be forced
          anytime by calling :func:`get_configuration`.
          """

          self.filename = kwargs.pop('filename', CONFIGNAME)

          args = tuple([self] + list(args))
          SafeConfigParser.__init__(*args, **kwargs)  # old style class ... grmbl
          # defaults
          self.set('DEFAULT', 'configdir', defaults['configdir'])
          self.set('DEFAULT', 'qscriptdir',
                  os.path.join("%(configdir)s", os.path.basename(defaults['qscriptdir'])))
          self.set('DEFAULT', 'templatesdir',
                  os.path.join("%(configdir)s", os.path.basename(defaults['templatesdir'])))
          self.set('DEFAULT', 'managerdir',
                  os.path.join("%(configdir)s", os.path.basename(defaults['managerdir'])))
          self.add_section('Gromacs')
          self.set("Gromacs", "tools", "pdb2gmx editconf grompp genbox genion mdrun trjcat trjconv")
          self.set("Gromacs", "extra", "")
          self.set("Gromacs", "groups", "tools")
          self.add_section('Logging')
          self.set('Logging', 'logfilename', defaults['logfilename'])
          self.set('Logging', 'loglevel_console', defaults['loglevel_console'])
          self.set('Logging', 'loglevel_file', defaults['loglevel_file'])

          # bundled defaults (should be ok to use get_template())
          default_cfg = get_template(self.cfg_template)
          self.readfp(open(default_cfg))

          # defaults are overriden by existing user global cfg file
          self.read([self.filename])

     @property
     def configuration(self):
          """Dict of variables that we make available as globals in the module.

          Can be used as ::

             globals().update(GMXConfigParser.configuration)        # update configdir, templatesdir ...
          """
          configuration = {
               'configfilename': self.filename,
               'logfilename': self.getpath('Logging', 'logfilename'),
               'loglevel_console': self.getLogLevel('Logging', 'loglevel_console'),
               'loglevel_file': self.getLogLevel('Logging', 'loglevel_file'),
               'configdir': self.getpath('DEFAULT', 'configdir'),
               'qscriptdir': self.getpath('DEFAULT', 'qscriptdir'),
               'templatesdir': self.getpath('DEFAULT', 'templatesdir'),
               }
          configuration['path'] = [os.path.curdir,
                                   configuration['qscriptdir'],
                                   configuration['templatesdir']]
          return configuration

     def getpath(self, section, option):
          """Return option as an expanded path."""
          return os.path.expanduser(os.path.expandvars(self.get(section, option)))

     def getlist(self, section, option, separator=None, sort=True):
          """Return *option* as a a list of strings.

          Split on white space or *separator*. Sort for *sort* = ``True``.
          """
          l = self.get(section, option).split(separator)
          if sort:
               l.sort()
          return l

     def getLogLevel(self, section, option):
          """Return the textual representation of logging level 'option' or the number.

          Note that option is always interpreted as an UPPERCASE string
          and hence integer log levels will not be recognized.

          .. SeeAlso: :mod:`logging` and :func:`logging.getLevelName`
          """
          return logging.getLevelName(self.get(section, option).upper())

def get_configuration(filename=CONFIGNAME):
    """Reads and parses the configuration file.

    Default values are loaded and then replaced with the values from
    ``~/.gromacswrapper.cfg`` if that file exists. The global
    configuration instance :data:`gromacswrapper.config.cfg` is updated
    as are a number of global variables such as :data:`configdir`,
    :data:`qscriptdir`, :data:`templatesdir`, :data:`logfilename`, ...

    Normally, the configuration is only loaded when the :mod:`gromacs`
    package is imported but a re-reading of the configuration can be forced
    anytime by calling :func:`get_configuration`.

    :Returns: a dict with all updated global configuration variables
    """
    global cfg, configuration    # very iffy --- most of the whole config mod should a class

    #: :data:`cfg` is the instance of :class:`GMXConfigParser` that makes all
    #: global configuration data accessible
    cfg = GMXConfigParser(filename=filename)   # update module-level cfg
    globals().update(cfg.configuration)        # update configdir, templatesdir ...
    configuration = cfg.configuration          # update module-level configuration
    return cfg

#: :data:`cfg` is the instance of :class:`GMXConfigParser` that makes all
#: global configuration data accessible
cfg = GMXConfigParser()
globals().update(cfg.configuration)        # update configdir, templatesdir ...

#: Dict containing important configuration variables, populated by
#: :func:`get_configuration` (mainly a shortcut; use :data:`cfg` in most cases).
configuration = cfg.configuration

def setup(filename=CONFIGNAME):
     """Prepare a default GromacsWrapper global environment.

     1) Create the global config file.
     2) Create the directories in which the user can store template and config files.

     This function can be run repeatedly without harm.
     """
     # setup() must be separate and NOT run automatically when config
     # is loaded so that easy_install installations work
     # (otherwise we get a sandbox violation)
     # populate cfg with defaults (or existing data)
     get_configuration()
     if not os.path.exists(filename):
          with open(filename, 'w') as configfile:
               cfg.write(configfile)  # write the default file so that user can edit
               msg = "NOTE: GromacsWrapper created the configuration file \n\t%r\n" \
                     "      for you. Edit the file to customize the package." % filename
               print msg

     # directories
     for d in config_directories:
          utilities.mkdir_p(d)


def check_setup():
     """Check if templates directories are setup and issue a warning and help.

     Returns ``True`` if all files and directories are found and
     ``False`` otherwise.

     Setting the environment variable
     :envvar:`GROMACSWRAPPER_SUPPRESS_SETUP_CHECK` to 'true' ('yes'
     and '1' also work) silence this function and make it always return ``True``.

     .. versionchanged:: 0.3.1
        Uses :envvar:`GROMACSWRAPPER_SUPPRESS_SETUP_CHECK` to suppress output
        (useful for scripts run on a server)
     """
     if os.environ.get("GROMACSWRAPPER_SUPPRESS_SETUP_CHECK", "false").lower() in ("1", "true", "yes"):
         return True

     is_complete = True
     if not os.path.exists(CONFIGNAME):
          is_complete = False
          print "NOTE: The global configuration file %(CONFIGNAME)r is missing." % globals()
          print "      You can create a basic configuration file from within python:"
          print "      >>> import gromacs"
          print "      >>> gromacs.config.setup()"

     missing = [d for d in config_directories if not os.path.exists(d)]
     if len(missing) > 0:
          is_complete = False
          print "NOTE: Some configuration directories are not set up yet"
          print "      %r" % missing
          print "      You can create them with the command from within Python"
          print "      >>> import gromacs"
          print "      >>> gromacs.config.setup()"
     return is_complete

check_setup()



# Gromacs tools
# -------------

#: Python list of all tool file names. Filled from values in the tool
#: groups in the configuration file.
load_tools = []
for g in cfg.getlist('Gromacs', 'groups', sort=False):
     load_tools.extend(cfg.getlist('Gromacs', g))
del g

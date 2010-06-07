# GromacsWrapper config.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
:mod:`gromacs.config` -- Configuration for GromacsWrapper
==========================================================

The config module provides configurable options for the whole package;
eventually it might grow into a sophisticated configuration system such as
matplotlib's rc system but right now it mostly serves to define which gromacs
tools and other scripts are offered in the package and where template files are
located. If the user wants to change anything they will still have to do it
here in source until a better mechanism with rc files has been implemented.

User-supplied templates are stored under
:data:`gromacs.config.configdir`. Eventually this will also contain the
configuration options currently hard-coded in :mod:`gromacs.config`.

.. autodata:: configdir
.. autodata:: path

The user should execute :func:`gromacs.config.setup` at least once to
prepare the user configurable area in their home directory::

  import gromacs
  gromacs.config.setup()


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
.. autodata:: load_scripts

:data:`load_tools` is populated by listing ``gmx_*`` tool group variables in
:data:`gmx_tool_groups`. 

.. autodata:: gmx_tool_groups

The tool groups variables are strings that contain white-space separated file
names of Gromacs tools. These lists determine which tools are made available as
classes in :mod:`gromacs.tools`.

.. autodata:: gmx_tools
.. autodata:: gmx_extra_tools


Location of template files
--------------------------

*Template variables* list files in the package that can be used as
templates such as run input files. Because the package can be a zipped
egg we actually have to unwrap these files at this stage but this is
completely transparent to the user.

.. autodata:: qscriptdir
.. autodata:: templatesdir
.. autodata:: templates
.. autodata:: qscript_template
.. autofunction:: setup

Accessing configuration data
----------------------------

The following functions can be used to access configuration data. Note that
files are searched first with their full filename, then in all directories
listed in :data:`gromacs.config.path`, and finally within the package itself.

.. autofunction:: get_template
.. autofunction:: get_templates

"""

import os
from pkg_resources import resource_filename, resource_listdir

import utilities

# Logging
# -------
import logging

#: File name for the log file; all gromacs command and many utility functions (e.g. in
#: :mod:`gromacs.cbook` and :mod:`gromacs.setup`) append messages there. Warnings and
#: errors are also recorded here. The default is *gromacs.log*.
logfilename = "gromacs.log"

#: The default loglevel that is still printed to the console.
loglevel_console = logging.INFO

#: The default loglevel that is still written to the :data:`logfilename`.
loglevel_file = logging.DEBUG


# User-accessible configuration
# -----------------------------

#: Directory to store user templates and rc files.
#: The default value is ``~/.gromacswrapper``.
configdir = os.path.expanduser(os.path.join("~",".gromacswrapper"))

#: Directory to store user supplied queuing system scripts.
#: The default value is ``~/.gromacswrapper/qscripts``.
qscriptdir = os.path.join(configdir, 'qscripts')

#: Directory to store user supplied template files such as mdp files.
#: The default value is ``~/.gromacswrapper/templates``.
templatesdir = os.path.join(configdir, 'templates')

def setup():
     """Create the directories in which the user can store template and config files.

     This function can be run repeatedly without harm.
     """
     # setup() must be separate and NOT run automatically when config
     # is loaded so that easy_install installations work
     # (otherwise we get a sandbox violation)
     utilities.mkdir_p(configdir)
     utilities.mkdir_p(qscriptdir)
     utilities.mkdir_p(templatesdir)

def check_setup():
     """Check if templates directories are setup and issue a warning and help."""
     missing = [d for d in (configdir, qscriptdir, templatesdir)
                if not os.path.exists(d)]
     if len(missing) > 0:
          print "NOTE: Some configuration directories are not set up yet"
          print "      %r" % missing
          print "      You can create them with the command"
          print "      >>> gromacs.config.setup()"
     return len(missing) == 0
check_setup()
               

#: Search path for user queuing scripts and templates. The internal package-supplied
#: templates are always searched last via :func:`gromacs.config.get_templates`. 
#: Modify :data:`gromacs.config.path` directly in order to customize the template 
#: and qscript searching. By default it has the value ``['.', qscriptdir,
#: templatesdir]``. 
#: (Note that it is not a good idea to have template files and qscripts with the
#: same name as they are both searched on the same path.)
path = [os.path.curdir, qscriptdir, templatesdir]

# Gromacs tools
# -------------

#: List of the variables in gromacs.tools that should be loaded. Possible values:
#: *gmx_tools*, *gmx_extra_tools*. Right now these are variable names in
#: :mod:`gromacs.config`, referencing :data:`gromacs.config.gmx_tools` and
#: :data:`gromacs.config.gmx_extra_tools`.
gmx_tool_groups = ['gmx_tools', 'gmx_extra_tools' ]

#: Contains the file names of all Gromacs tools for which classes are generated.
#: Editing this list has only an effect when the package is reloaded.  If you want
#: additional tools then add the, to the source (``config.py``) or derive new
#: classes manually from :class:`gromacs.core.GromacsCommand`.  (Eventually, this
#: functionality will be in a per-user configurable file.)  The current list was
#: generated from Gromacs 4.0.99 (git).
#: Removed (because of various issues)
#:   - g_kinetics
gmx_tools = """\
anadock      g_current     g_helix        g_rama     g_traj     mdrun_d
             g_density     g_helixorient  g_rdf      g_vanhove  mk_angndx
do_dssp      g_densmap                    g_rms      g_velacc   pdb2gmx
editconf     g_dielectric  g_lie                     g_wham     protonate
eneconv      g_dih         g_mdmat        g_rmsf     genbox     sigeps
g_anaeig     g_dipoles     g_mindist      g_rotacf   genconf    tpbconv
g_analyze    g_disre       g_morph        g_saltbr   genion     trjcat
g_angle      g_dist        g_msd          g_sas      genrestr   trjconv
g_bond       g_dyndom      g_nmeig        g_sdf      gmxcheck   trjorder
g_bundle     g_enemat      g_nmens        g_sgangle  gmxdump    wheel
g_chi        g_energy      g_nmtraj       g_sham     grompp     x2top
g_cluster    g_filter      g_order        g_sorient  luck       
g_clustsize  g_gyrate      g_polystat     g_spatial  make_edi   xpm2ps
g_confrms    g_h2order     g_potential    g_spol     make_ndx
g_covar      g_hbond       g_principal    g_tcaf     mdrun
"""

gmx_tools_402 = """\
anadock      g_current     g_helix        g_rama     g_traj     mdrun_d
demux.pl     g_density     g_helixorient  g_rdf      g_vanhove  mk_angndx
do_dssp      g_densmap     g_kinetics     g_rms      g_velacc   pdb2gmx
editconf     g_dielectric  g_lie          g_rmsdist  g_wham     protonate
eneconv      g_dih         g_mdmat        g_rmsf     genbox     sigeps
g_anaeig     g_dipoles     g_mindist      g_rotacf   genconf    tpbconv
g_analyze    g_disre       g_morph        g_saltbr   genion     trjcat
g_angle      g_dist        g_msd          g_sas      genrestr   trjconv
g_bond       g_dyndom      g_nmeig        g_sdf      gmxcheck   trjorder
g_bundle     g_enemat      g_nmens        g_sgangle  gmxdump    wheel
g_chi        g_energy      g_nmtraj       g_sham     grompp     x2top
g_cluster    g_filter      g_order        g_sorient  luck       xplor2gmx.pl
g_clustsize  g_gyrate      g_polystat     g_spatial  make_edi   xpm2ps
g_confrms    g_h2order     g_potential    g_spol     make_ndx
g_covar      g_hbond       g_principal    g_tcaf     mdrun
"""


#: Additional gromacs tools (add *gmx_extra_tools* to
#: :data:`gromacs.config.gmx_tool_groups` to enable them, provided the binaries
#: have been provided on the :envvar:`PATH`).
gmx_extra_tools = """\
g_count      g_flux        g_zcoord
g_ri3Dc      a_ri3Dc       a_gridcalc
"""

#: Python list of all tool file names. Automatically filled from :data:`gmx_tools`
#: and :data:`gmx_extra_tools`, depending on the values in
#: :data:`gmx_tool_groups`.
load_tools = []

for g in gmx_tool_groups:
     load_tools.extend(globals()[g].split())
del g

# Adding additional 3rd party scripts
# -----------------------------------

# HACK-Alert: we're temporarily "installing" scripts from the egg  ... this is a hack:
# Must extract because it is part of a zipped python egg;
# see http://peak.telecommunity.com/DevCenter/PythonEggs#accessing-package-resources
GridMAT_MD = resource_filename(__name__,'external/GridMAT-MD_v1.0.2/GridMAT-MD.pl')
os.chmod(GridMAT_MD, 0755)


#: 3rd party analysis scripts and tools; this is a list of triplets of
#:
#:   (*script name/path*, *command name*, *doc string*)
#:
#: (See the source code for examples.) 
load_scripts = [
    (GridMAT_MD,
     'GridMAT_MD',
     """GridMAT-MD: A Grid-based Membrane Analysis Tool for use with Molecular Dynamics.

*This* ``GridMAT-MD`` is a patched version of the original
``GridMAT-MD.pl`` v1.0.2, written by WJ Allen, JA Lemkul and DR
Bevan. The original version is available from the `GridMAT-MD`_ home
page,
   
.. _`GridMAT-MD`: http://www.bevanlab.biochem.vt.edu/GridMAT-MD/index.html

Please cite 

  W. J. Allen, J. A. Lemkul, and D. R. Bevan. (2009) "GridMAT-MD: A
  Grid-based Membrane Analysis Tool for Use With Molecular Dynamics."
  J. Comput. Chem. 30 (12): 1952-1958.

when using this programme.

Usage:

.. class:: GridMAT_MD(config[,structure])

:Arguments:
   - *config* : See the original documentation for a description for the
     configuration file.
   - *structure* : A gro or pdb file that overrides the value for
     *bilayer* in the configuration file.

"""),
    ]


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


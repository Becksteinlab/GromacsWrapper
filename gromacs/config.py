# $Id$
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

.. autodata:: templates
.. autodata:: sge_template


Functions
---------

The following functions can be used to access configuration data.

.. autofunction:: get_template

"""

import os
from pkg_resources import resource_filename


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
loglevel_file = logging.INFO


# Gromacs tools
# -------------

#: List of the variables in gromacs.tools that should be loaded. Possible values:
#: *gmx_tools*, *gmx_extra_tools*. Right now these are variable names in
#: :mod:`gromacs.config`, referencing data:`gromacs.config.gmx_tools` and
#: data:`gromacs.config.gmx_extra_tools`.
gmx_tool_groups = ['gmx_tools', ]

#: Contains the file names of all Gromacs tools for which classes are generated.
#: Editing this list has only an effect when the package is reloaded.  If you want
#: additional tools then add the, to the source (``config.py``) or derive new
#: classes manually from :class:`gromacs.core.GromacsCommand`.  (Eventually, this
#: functionality will be in a per-user configurable file.)  The current list was
#: generated from Gromacs 4.0.2.
gmx_tools = """\
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

# XXX: install script via setup.py ... this is a hack:
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

# TODO: This is becoming unwieldy: should be more abstract/automatic.
#       Maybe simply use filename as key? This would require some cleaning up
#       in gromacs.setup but would not be too bad.
templates = {
    'em_mdp': resource_filename(__name__, 'templates/em.mdp'),
    'md_G43a1_mdp': resource_filename(__name__, 'templates/md_G43a1.mdp'),
    'md_OPLSAA_mdp': resource_filename(__name__, 'templates/md_OPLSAA.mdp'),
    'deathspud_sge': resource_filename(__name__, 'templates/deathspud.sge'),
    'neuron_sge': resource_filename(__name__, 'templates/neuron.sge'),
    'hector_pbs': resource_filename(__name__, 'templates/hector.pbs'),
    'hpcx_ll': resource_filename(__name__, 'templates/hpcx.ll'),    
    }
"""Templates have to be extracted from the egg because they are used
by external code. All template filenames are stored in
:data:`gromacs.config.templates`.

**Gromacs mdp templates**

   These are supplied as examples and there is *NO GUARANTEE THAT THEY
   PRODUCE SENSIBLE OUTPUT* --- check for yourself!  Note that only
   existing parameter names can be modified with
   :func:`gromacs.cbook.edit_mdp` at the moment; if in doubt add the
   parameter with its gromacs default value (or empty values) and
   modify later with :func:`~gromacs.cbook.edit_mdp`.

**SGE templates**

   (This is a misnomer --- there are also scripts for PBS and
   LoadLeveller but we call them all SGE scripts... apologies).

   The sge scripts are highly specific and you will need to add your
   own.  Templates should be sh-scripts and can contain the following
   patterns; these are either shell variable assignments or batch
   submission system commands. The table shows SGE commands but PBS
   and LoadLeveller have similar constructs; e.g. PBS commands start
   with ``#PBS`` and LoadLeveller uses ``#@`` with its own comman
   keywords):


   ===============  ===========  ================ ================= =====================================
   command          default      replacement      description       regex
   ===============  ===========  ================ ================= =====================================
   #$ -N            GMX_MD       *sgename*        job name          /^#.*(-N|job_name)/
   #$ -l walltime=  00:20:00     *walltime*       max run time      /^#.*(-l walltime|wall_clock_limit)/
   #$ -A            BUDGET       *budget*         account           /^#.*(-A|account_no)/
   DEFFNM=          md           *deffnm*         default gmx name  /^DEFFNM=/
   WALL_HOURS=      0.33         *walltime* h     mdrun's -maxh     /^WALL_HOURS=/
   ===============  ===========  ================ ================= =====================================

   These lines should not have any white space at the beginning. The
   regular expression pattern is used to find the lines for the
   replacement and the default values are replaced.

"""

#: The default template for SGE/PBS run scripts.
sge_template = templates['neuron_sge']


# Functions to access configuration data
# --------------------------------------

def get_template(t):
    """Find template file *t* and return its real path.

    *t* can be

    1. a relative or absolute path,
    2. a filename in the package template directory (defined in the template dictionary
       :data:`gromacs.config.templates`) or
    3. a key into :data:`~gromacs.config.templates`.

    The first match (in this order) is returned.

    :Arguments: *t* : template file or key
    :Returns:   os.path.realpath(*t*)
    :Raises:    :exc:`ValueError` if no file can be located.
       
    """
    if os.path.exists(t):           # 1) Is it an accessible file?
        pass
    else:                           # 2) check the packaged template files
        _t = os.path.basename(t)
        _t_found = False
        for p in templates.values():
            if _t == os.path.basename(p):
                t = p
                _t_found = True     # NOTE: in principle this could match multiple
                break               #       times if more than one template dir existed.
        if not _t_found:            # 3) try it as a key into templates
            try:
                t = templates[t]
            except KeyError:
                pass
            else:
                _t_found = True
        if not _t_found:            # 4) nothing else to try... or introduce a PATH?
            raise ValueError("Failed to locate the template file %(t)r." % vars())
    return os.path.realpath(t)


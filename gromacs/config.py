# $Id$
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
:mod:`gromacs.config` -- Configuration for GromacsWrapper
==========================================================

The config module is supposed to provide configurable options for the
whole package; eventually it might grow into a sophisticated
configuration system such as matplotlib's rc system but right now it
only serves to define which gromacs tools and other scripts are
offered in the package. If the user wants to change anything they will
still have to do it here (and in :mod:`gromacs.tools`) until a better
mechanism with rc files has been implemented.

Data
----

``load_*`` variables are lists that contain instructions to other
parts of the code which packages and scripts should be wrapped.

.. autodata:: load_tools
.. autodata:: load_scripts

Template variables list files in the package that can be used as
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


#: preliminary: list the variables in gromacs.tools that should be
#: loaded. Possible values: *gmx_tools*, *gmx_extra_tools*.
load_tools = ['gmx_tools', ]

# XXX: install script via setup.py ... this is a hack:
# Must extract because it is part of a zipped python egg;
# see http://peak.telecommunity.com/DevCenter/PythonEggs#accessing-package-resources
GridMAT_MD = resource_filename(__name__,'external/GridMAT-MD_v1.0.2/GridMAT-MD.pl')
os.chmod(GridMAT_MD, 0755)


#: 3rd party analysis scripts and tools; triplets of 
#:   (script name/path, command name, doc string)
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


templates = {
    'em_mdp': resource_filename(__name__, 'templates/em.mdp'),
    'md_G43a1_mdp': resource_filename(__name__, 'templates/md_G43a1.mdp'),
    'md_OPLSAA_mdp': resource_filename(__name__, 'templates/md_OPLSAA.mdp'),
    'deathspud_sge': resource_filename(__name__, 'templates/deathspud.sge'),
    'neuron_sge': resource_filename(__name__, 'templates/neuron.sge'),
    }
"""Templates have to be extracted from the egg because they are used
by external code. All template filenames are stored in
:data:`gromacs.config.templates`.

**Gromacs mdp templates**

   These are supplied as examples and there is *NO GUARANTEE THAT THEY
   PRODUCE SENSIBLE OUTPUT* --- check for yourself!  Note that only
   existing parameter names can be modified with
   :func:`gromacs.cbook.edit_mdp` at the moment; if in doubt add the the
   parameter with its gromacs default value (or empty values) and
   modify later with :func:`~gromacs.cbook.edit_mdp`.

**SGE templates**

   The sge scripts are highly specific and you will need to add your own.
   Templates should be sh-scripts and contain the following lines::
      #$ -N GMX_MD     'GMX_MD' is replaced by kw sgename
      DEFFNM=md        'md' is replaced by kw deffnm
"""

#: The default template for SGE run scripts.
sge_template = templates['neuron_sge']


def get_template(t):
    """Find template file *t* and return its real path.

    *t* can be a relative or absolute path, a filename in the package template
    directory (defined in the template dictionary
    :data:`gromacs.config.templates`) or a key into
    :data:`~gromacs.config.templates`. The first match (in this order) is
    returned.

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


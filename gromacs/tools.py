# $Id$
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
gromacs.tools -- Gromacs commands classes
=========================================

A Gromacs command class can be thought of as a factory function that
produces an instance of a gromacs command with initial default values.

By convention, a class has the capitalized name of the corresponding Gromacs
tool; dots are replaced by underscores to make it a valid python identifier.

At the moment the gromacs tools are hard coded in ``gromacs.tools.gmx_tools``;
the list was generated from Gromacs 4.0.2.

.. autodata:: gmx_tools
.. autodata:: gmx_extra_tools

Example
-------

In this example we create two instances of the ``Trjconv`` command::

  import gromacs.tools as tools

  trjconv = tools.Trjconv()
  trjconv_compact = tools.Trjconv(ur='compact', center=True, boxcenter='tric', pbc='mol',
                                  input=('protein','system'),
                                  doc="Returns a compact representation of the system centered on the protein")

The first one ``trjconv`` behaves as the standard commandline tool but the
second one, ``trjconv_compact``, will by default create a compact
representation of the input data by taking into account the shape of the unit
cell. Of course, the same effect can be obtained by providing the corresponding
arguments to ``trjconv`` but by naming the more specific command differently
one can easily build up a library of small tools that will solve a specifi,
repeatedly encountered problem reliably. This is particularly helpful when doin
interactive work.
"""

__docformat__ = "restructuredtext en"

from core import GromacsCommand

# Auto-generate classes such as:
# class g_dist(GromacsCommand):
#     command_name = 'g_dist'

#: Contains the file names of all Gromacs tools for which classes are generated.
#: Changing this list does *not* add additional classes. Either change the source
#: or derive new classes manually from ``gromacs.core.GromacsCommand``.
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

#: Additional gromacs tools (not added at the moment).
gmx_extra_tools = """\
g_count      g_flux
g_ri3Dc      a_ri3Dc       a_gridcalc
"""

registry = {}

for name in gmx_tools.split():
    # make names valid python identifiers and convention: class names are capitalized
    clsname = name.replace('.','_').capitalize()  
    cls = type(clsname, (GromacsCommand,), {'command_name':name})
    registry[clsname] = cls      # registry keeps track of all classes

locals().update(registry)        # add classes to module's scope

del name, cls, clsname

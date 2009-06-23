# $Id$
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""Gromacs commands classes.

A class can be thought of as a factory function that produces an
instance of a gromacs command with initial default values.

A class has the capitalized name of the corresponding Gromacs tool;
dots are replaced by underscores to make it a valid python identifier.
"""
__docformat__ = "restructuredtext en"

from core import GromacsCommand

# Auto-generate classes such as:
# class g_dist(GromacsCommand):
#     command_name = 'g_dist'

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

# not added at the moment
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

del name, cls, clsname, gmx_tools

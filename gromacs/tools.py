# $Id$
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
:mod:`gromacs.tools` -- Gromacs commands classes
================================================

A Gromacs command class can be thought of as a factory function that
produces an instance of a gromacs command with initial default values.

By convention, a class has the capitalized name of the corresponding Gromacs
tool; dots are replaced by underscores to make it a valid python identifier.

At the moment the gromacs tools are hard coded in :data:`gromacs.tools.gmx_tools`;
the list was generated from Gromacs 4.0.2.

.. autodata:: gmx_tools
.. autodata:: gmx_extra_tools

It is also possible to extend the basic commands and patch in additional
functionality. For example, the :class:`GromacsCommandMultiIndex` class makes a
command accept multiple index files and concatenates them on the fly; the
behaviour mimics Gromacs' "multi-file" input that has not yet been enabled for
all tools.

.. autoclass:: GromacsCommandMultiIndex
   :members: __init__, run, _fake_multi_index, __del__

Example
-------

In this example we create two instances of the :class:`gromacs.tools.Trjconv` command (which
runs the Gromacs ``trjconv`` command)::

  import gromacs.tools as tools

  trjconv = tools.Trjconv()
  trjconv_compact = tools.Trjconv(ur='compact', center=True, boxcenter='tric', pbc='mol',
                                  input=('protein','system'),
                                  doc="Returns a compact representation of the system centered on the protein")

The first one, ``trjconv``, behaves as the standard commandline tool but the
second one, ``trjconv_compact``, will by default create a compact
representation of the input data by taking into account the shape of the unit
cell. Of course, the same effect can be obtained by providing the corresponding
arguments to ``trjconv`` but by naming the more specific command differently
one can easily build up a library of small tools that will solve a specifi,
repeatedly encountered problem reliably. This is particularly helpful when doing
interactive work.
"""

__docformat__ = "restructuredtext en"

import tempfile

from core import GromacsCommand
import utilities

#: Contains the file names of all Gromacs tools for which classes are generated.
#: Changing this list does *not* add additional classes. Either change the source
#: or derive new classes manually from :class:`gromacs.core.GromacsCommand`.
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

#: This dict holds all generated classes.
registry = {}

# Auto-generate classes such as:
# class g_dist(GromacsCommand):
#     command_name = 'g_dist'

for name in gmx_tools.split():
    # make names valid python identifiers and convention: class names are capitalized
    clsname = name.replace('.','_').capitalize()  
    cls = type(clsname, (GromacsCommand,), {'command_name':name,
                                            '__doc__': "Gromacs tool %(name)r." % vars()})
    registry[clsname] = cls      # registry keeps track of all classes

globals().update(registry)        # add classes to module's scope

del name, cls, clsname

# modify/fix classes as necessary
# Note: 
# - check if class was defined in first place
# - replace class
# - update local context AND registry as done below

class GromacsCommandMultiIndex(GromacsCommand):
        def __init__(self, **kwargs):
            kwargs = self._fake_multi_ndx(**kwargs)
            super(GromacsCommandMultiIndex, self).__init__(**kwargs)

        def run(self,*args,**kwargs):
            kwargs = self._fake_multi_ndx(**kwargs)            
            return super(GromacsCommandMultiIndex, self).run(*args, **kwargs)

        def _fake_multi_ndx(self, **kwargs):
            """Combine multiple index file into a single one and return appropriate kwargs.

            Calling the method combines multiple index files into a a single
            temporary one so that Gromacs tools that do not (yet) support multi
            file input for index files can be used transparently as if they did.

            If a temporary index file is required then it is deleted once the
            object is destroyed.

            :Returns:
            The method returns the input keyword arguments with the necessary
            changes to use the temporary index files.

            :Keywords: Only the listed keywords have meaning for the method:
               n : filename or list of filenames
                  possibly multiple index files; *n* is replaced by the name of 
                  the temporary index file.
               s : filename
                  structure file (tpr, pdb, ...) or ``None``; if a structure file is 
                  supplied then the Gromacs default index groups are automatically added
                  to the temporrary indexs file.

            :Example: 
               Used in derived classes that replace the standard
               :meth:`run` (or :meth:`__init__`) methods with something like::

                  def run(self,*args,**kwargs):
                      kwargs = self._fake_multi_ndx(**kwargs)            
                      return super(G_mindist, self).run(*args, **kwargs)

                      """
            ndx = kwargs.get('n')
            if not (ndx is None or type(ndx) is str):
                if len(ndx) > 1:
                    # g_mindist cannot deal with multiple ndx files (at least 4.0.5)
                    # so we combine them in a temporary file; it is unlinked in __del__.
                    # self..multi_ndx stores file name for __del__
                    fd, self.multi_ndx = tempfile.mkstemp(suffix='.ndx', prefix='multi_')
                    make_ndx = Make_ndx(f=kwargs.get('s'), n=ndx)
                    rc,out,err = make_ndx(o=self.multi_ndx, input=['q'], 
                                          stdout=False, stderr=False)
                    self.orig_ndx = ndx
                    kwargs['n'] = self.multi_ndx
            return kwargs

        def __del__(self):
            """Clean up temporary multi-index files if they were used."""
            try:
                # self.multi_ndx <-- _fake_multi_index()
                utilities.unlink_gmx(self.multi_ndx)
            except (AttributeError, OSError):
                pass
            # XXX: type error --- can't use super in __del__?
            #super(GromacsCommandMultiIndex, self).__del__()



if 'G_mindist' in registry:
    del G_mindist, registry['G_mindist']

    # let G_mindist handle multiple ndx files
    class G_mindist(GromacsCommandMultiIndex):
        """Gromacs tool 'g_mindist' (with patch to handle multiple ndx files)."""
        command_name = 'g_mindist'

    registry['g_mindist'] = G_mindist


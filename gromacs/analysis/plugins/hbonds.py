# -*- encoding: utf-8 -*-
# Copyright (c) 2012 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
Hydrogen bond analysis
======================

Analysis of hydrogen bonds using :class:`gromacs.tools.G_hbond`.

Example
-------

Analysis of hydrogen bonds between a ligand (resname 5FH) and the
protein/water::

  import gromacs.analysis
  gromacs.start_logging()
  S = gromacs.analysis.Simulation(tpr="MD/md.tpr", xtc="MD/md.xtc", analysisdir="./analysis", \
      plugins=[gromacs.analysis.plugins.HBonds(group1='r 5FH', group2='"Protein" | "Water"')])
  S.set_plugin('HBonds')
  S.run('HBonds', force=True, smooth=-1)
  S.analyze()

* Output files will be produced under the analysisdir, ``./analysis/HBonds``,
  including the file ``hb_existence.txt`` which lists the percentage of the
  total simulation during which the hydrogen bond was observed.

* The selection of groups is specified via :class:`gromacs.tools.Make_ndx`
  syntax (in fact, :class:`gromacs.cbook.IndexBuilder` is used but currently
  the special command syntax of :class:`~gromacs.cbook.IndexBuilder` is not
  enabled here).


Plugin class
------------

.. autoclass:: HBonds
   :members: worker_class
   :undoc-members:

Worker class
------------

The worker class performs the analysis.

.. autoclass:: _HBonds
   :members:


"""
from __future__ import with_statement

__docformat__ = "restructuredtext en"

import os.path
import warnings

import numpy as np

import gromacs
import gromacs.tools
from gromacs.utilities import AttributeDict
from gromacs.analysis.core import Worker, Plugin

import logging
logger = logging.getLogger('gromacs.analysis.plugins.HBonds')

# Worker classes that are registered via Plugins (see below)
# ----------------------------------------------------------
# These must be defined before the plugins.

class _HBonds(Worker):
    """HBonds worker class."""

    def __init__(self, **kwargs):
        """Set up  HBonds analysis.

        Hydrogen bonds Donor-H...Acceptor are analyzed. Depending on the
        *N_is_acceptor* flag, nitrogens can be counted as donors *and* acceptors.

        :Arguments:
           *group1*
               :func:`gromacs.make_ndx` selection string, e.g. 'r 5FH' ["Protein"]
           *group2*
               :func:`gromacs.make_ndx` selection string, e.g. '"Protein"'
               (note the double quotes) ["Protein"]
           *dmax*
               maximum length of an hydrogen bond (donor-acceptor), can be changed
               with *use_donor_acceptor_distance* [0.35 nm]
           *amax*
               maximum deviation of the acceptor-donor-hydrogen angle from
               linearity (which is 0º) [30º]
           *use_donor_acceptor_distance*
               ``True``: *dmax* refers to the donor-acceptor distance;
               ``False``: *dmax* is the hydrogen-acceptor max distance. [``True``]
           *N_is_acceptor*
               Count nitrogens as donors *and* acceptors (generally only
               important for certain nitrogens in heterocycles as in histidine)
               ? [``False``]

        """
        # specific arguments: take them before calling the super class that
        # does not know what to do with them
        group1 = kwargs.pop('group1', '"Protein"')
        group2 = kwargs.pop('group2', '"Protein"')
        gmxargs = {'r': kwargs.pop('dmax', 0.35),
                   'a': kwargs.pop('amax', 30.0),
                   'da': kwargs.pop('use_donor_acceptor_distance', True),
                   'nitacc': kwargs.pop('N_is_acceptor', False),
                   }

        # super class init: do this before doing anything else
        # (also sets up self.parameters and self.results)
        super(_HBonds, self).__init__(**kwargs)

        self.parameters.groups = {'group1': '@'+group1,
                                  'group2': '@'+group2,
                                  }  #  '@' for IndexBuilder raw command escape!!

        # if fit not good then g_hbond returns non-zero return code: just warn;
        # also pre-set the defaults
        self.g_hbond = gromacs.tools.G_hbond(failure="warn", **gmxargs)
        self.default_g_hbond_args = gmxargs

        # self.simulation might have been set by the super class
        # already; just leave this snippet at the end. Do all
        # initialization that requires the simulation class in the
        # _register_hook() method.
        if not self.simulation is None:
            self._register_hook()

    def _register_hook(self, **kwargs):
        """Run when registering; requires simulation."""

        super(_HBonds, self)._register_hook(**kwargs)
        assert not self.simulation is None

        self.parameters.filenames = {
            'ndx': self.plugindir('hb_groups.ndx'),          # filename of the index file
            'hbm': self.plugindir('hb.xpm'),
            'log': self.plugindir('hbond.log'),
            'num': self.plugindir('hbnum.xvg'),
            'dist': self.plugindir('hbdist.xvg'),
            'ang': self.plugindir('hbang.xvg'),
            'hbn': self.plugindir('hbond.ndx'),
            'ac':  self.plugindir('hbac.xvg'),
            'life': self.plugindir('hblife.xvg'),
            'existence': self.plugindir('hb_existence.txt'),
            }
        self.parameters.figname = self.figdir('hbonds')  # not used yet

    def create_ndx(self, ndx=None, **kwargs):
        """Create index file *ndx*.

        Uses *group1* and *group2* as selection commands in :func:`gromacs.make_ndx`.
        """
        if ndx is None:
            ndx = self.parameters.filenames['ndx']
        g = self.parameters.groups
        names = ['group1', 'group2']
        I = gromacs.cbook.IndexBuilder(struct=self.simulation.tpr, out_ndx=ndx,
                                       selections=[g[name] for name in names],
                                       names=names)
        I.write()
        return names

    def get_temperature_from_tpr(self, default=298.15):
        """Hack to get the temperature from the tpr file."""
        temperature = default
        rc,out,err = gromacs.gmxdump(s=self.simulation.tpr, stdout=False, stderr=False)
        for line in out.splitlines():
            if line.strip().startswith('ref_t'):
                temperature = np.mean(map(float, line[10:].split()))
                break
        return temperature

    def run(self, force=False, **gmxargs):
        """Run g_hbonds.

        Processes the trajectory and writes the data files. Can take a
        long time if the solvent is involved.
        """
        F = self.parameters.filenames

        if not self.check_file_exists(F['ndx'], resolve='warning') or force:
            logger.info("Creating the index groups with commands %r", self.parameters.groups)
            self.create_ndx(F['ndx'])

        if not self.check_file_exists(F['num'], resolve='warning') or force:
            temperature = self.get_temperature_from_tpr()
            gmxargs.setdefault('temp', temperature)
            gmxargs.setdefault('smooth', -1)
            logger.info("Analyzing HBonds (T=%g K from tpr)...", temperature)
            self.g_hbond(s=self.simulation.tpr, f=self.simulation.xtc, n=F['ndx'],
                         num=F['num'], g=F['log'], hbn=F['hbn'], hbm=F['hbm'], ac=F['ac'],
                         dist=F['dist'], ang=F['ang'], life=F['life'],
                         input=['group1', 'group2'],
                         **gmxargs)

    def analyze(self,**kwargs):
        """Analyze hydrogen bond output.

        * hydrogen bond existence (existence)
        * total number of hydrogen bonds (num)
        * (others can be added easily)

        :Returns:  a dictionary of the results and also sets ``self.results``.
        """
        from gromacs.formats import XPM, XVG

        results = AttributeDict()
        results['num'] = XVG(self.parameters.filenames['num'])
        results['matrix'] = hbm = XPM(self.parameters.filenames['hbm'], reverse=True)

        hb_fraction = hbm.array.mean(axis=0)
        desc = [line.strip() for line in
                open(self.parameters.filenames['log']) if not line.startswith('#')]
        results['existence'] = zip(desc, hb_fraction)

        with open(self.parameters.filenames['existence'], "w") as out:
            logger.info("Hydrogen bond existence analysis (results['existence'] and %(existence)r)",
                        self.parameters.filenames)
            for name,frac in results['existence']:
                logger.info("hb_existence: %-40s %4.1f%%", name, 100*frac)
                out.write("%-40s %4.1f%%\n" % (name, 100*frac))

        self.results = results
        return results

    def plot(self, **kwargs):
        """Plot all results, labelled by the result keys.

        :Keywords:
           figure
               - ``True``: save figures in the given formats
               - "name.ext": save figure under this filename (``ext`` -> format)
               - ``False``: only show on screen
           formats : sequence
               sequence of all formats that should be saved [('png', 'pdf')]
           plotargs
               keyword arguments for pylab.plot()

        TODO: plot individually, not all in one figure
        """

        import pylab
        figure = kwargs.pop('figure', False)
        extensions = kwargs.pop('formats', ('pdf','png'))
        for name,result in self.results.items():
            kwargs['label'] = name
            try:
                result.plot(**kwargs)      # This requires result classes with a plot() method!!
            except AttributeError:
                warnings.warn("Sorry, plotting of result %(name)r is not implemented" % vars(),
                              category=UserWarning)
        pylab.legend(loc='best')
        #if figure is True:
        #    for ext in extensions:
        #        self.savefig(ext=ext)
        #elif figure:
        #    self.savefig(filename=figure)




# Public classes that register the worker classes
#------------------------------------------------

class HBonds(Plugin):
    """*HBonds* plugin.

    Analysis of hydrogen bonds between two groups.

    Hydrogen bonds Donor-H...Acceptor are analyzed. Depending on the
    *N_is_acceptor* flag, nitrogens can be counted as donors *and* acceptors.


    .. class:: HBonds([group1, group2[, dmax, amax, ..., name[, simulation]]])

    :Arguments:
       *group1*
          :func:`gromacs.make_ndx` selection string, e.g. 'r 5FH' ["Protein"]
       *group2*
          :func:`gromacs.make_ndx` selection string, e.g. '"Protein"'
          (note the double quotes) ["Protein"]
       *dmax*
           maximum length of an hydrogen bond (donor-acceptor), can be changed
           with *use_donor_acceptor_distance* [0.35 nm]
       *amax*
           maximum deviation of the acceptor-donor-hydrogen angle from
           linearity (which is 0º) [30º]
       *use_donor_acceptor_distance*
           ``True``: *dmax* refers to the donor-acceptor distance;
           ``False``: *dmax* is the hydrogen-acceptor max distance. [``True``]
       *N_is_acceptor*
           Count nitrogens as donors *and* acceptors? [``True``]
        *name* : string
            plugin name (used to access it)
        *simulation* : instance
            The :class:`gromacs.analysis.Simulation` instance that owns the plugin.

    .. versionadded:: 0.3.1
    """
    worker_class = _HBonds



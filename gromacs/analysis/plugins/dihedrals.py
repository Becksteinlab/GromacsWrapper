# $Id$
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
Dihedrals plugin
================

Analyze one or more dihedral angles.

Plugin class
------------

.. autoclass:: Dihedrals
   :members: worker_class
   :undoc-members:

Worker class
------------

The worker class performs the analysis.

.. autoclass:: _Dihedrals
   :members:


"""
from __future__ import with_statement

__docformat__ = "restructuredtext en"

import os.path
import warnings

import numpy

import gromacs
from gromacs import MissingDataError
from gromacs.utilities import AttributeDict
from gromacs.analysis.core import Worker, Plugin
from gromacs.formats import XVG, NDX


# Worker classes that are registered via Plugins (see below)
# ----------------------------------------------------------
# These must be defined before the plugins.

class _Dihedrals(Worker):
    """Dihedrals worker class."""

    def __init__(self,**kwargs):
        """Set up dihedral analysis.

        :Arguments:
           *dihedrals*
             list of tuples; each tuple contains :class:`gromacs.cbook.IndexBuilder`
             atom selection commands.
           *labels*
             optional list of labels for the dihedrals. Must have as many entries as
             *dihedrals*.
        
        """
        # specific arguments: take them before calling the super class that
        # does not know what to do with them
        dihedralgroups = kwargs.pop('dihedrals', [])
        for dih in dihedralgroups:
            if len(dih) != 4 or type(dih) is str:
                raise ValueError("Each dihedral index group must contain exactly four atomnumbers, "
                                 "but this is not the case for %r." % dih)
        labels = kwargs.pop('labels', None)
        if not labels is None:
            if len(labels) != len(dihedralgroups):
                raise ValueError("Provide one label in labels for each dihedral in dihedrals.")

        # super class init: do this before doing anything else
        # (also sets up self.parameters and self.results)
        super(_Dihedrals, self).__init__(**kwargs)

        # process specific parameters now and set instance variables
        self.parameters.dihedrals = dihedralgroups
        self.parameters.labels = labels

        # self.simulation might have been set by the super class
        # already; just leave this snippet at the end. Do all
        # initialization that requires the simulation class in the
        # _register_hook() method.
        if not self.simulation is None:
            self._register_hook()

    def _register_hook(self, **kwargs):
        """Run when registering; requires simulation."""

        super(_Dihedrals, self)._register_hook(**kwargs)
        assert not self.simulation is None

        def P(*args):
            return self.plugindir(*args)

        self.parameters.ndx = P('dih.ndx')
        self.parameters.filenames = {
            'comb_distribution': P('dih_combdist.xvg'), # combined distribution from ALL dih!            
            'timeseries': P('dih_timeseries.xvg'),      # with -all, ie (phi, avg, dh1, dh2, ...)
            'transitions': P('dih_trans.xvg'),          # only works for multiplicity 3
            'transition_histogram': P('dih_transhisto.xvg'), 
            'acf': P('dih_acf.xvg'),
            'distributions' : P('dih_distributions.xvg'), # from analyze()
            'PMF': P('dih_pmf.xvg'),                      # from analyze()
            }
        self.parameters.figname = self.figdir('dihedrals_PMF')

        # just make index now... we have everything and it's really simple
        self.make_index()

    def make_index(self):
        """Make one index group *dihedrals* from the selections.

        Right now this requires **raw atom numbers** from the topology.
        """
        ndx = NDX()
        ndx['dihedrals'] = numpy.concatenate(self.parameters.dihedrals)
        ndx.write(filename=self.parameters.ndx, ncol=4)
        

    # override 'API' methods of base class
        
    def run(self, force=False, **gmxargs):
        """Collect dihedral data from trajectory with ``g_angle`` and save to data files.

        Note that the ``-all`` flag is always enabled and hence all time series
        can be found under the *timeseries* results key.
        """
        if not force and \
               self.check_file_exists(self.parameters.filenames['timeseries'], resolve='warn'):
            return

        angle_type = gmxargs.setdefault('type','dihedral')
        allowed_types = ("dihedral", "improper", "ryckaert-bellemans")
        if angle_type not in allowed_types:
            raise ValueError("The *type* can not be %r, only be one of\n\t%r\n" %
                             (angle_type, allowed_types))
        gmxargs['all'] = True   # required because we analyze the time series file
        gmxargs.setdefault('periodic', True)
        F = self.parameters.filenames
        gromacs.g_angle(f=self.simulation.xtc, n=self.parameters.ndx,
                        od=F['comb_distribution'], ov=F['timeseries'], ot=F['transitions'],
                        oc=F['acf'], oh=F['transition_histogram'], **gmxargs)


    def analyze(self, **kwargs):
        """Load results from disk into :attr:`_Dihedrals.results` and compute PMF.

        The PMF W(phi) in kT is computed from each dihedral
        probability distribution P(phi) as

           W(phi) = -kT ln P(phi)

        It is stored in :attr:`_Dihedrals.results` with the key *PMF*.

        :Keywords:
          *bins*
             bins for histograms (passed to numpy.histogram(new=True))
        
        :Returns: a dictionary of the results and also sets
                  :attr:`_Dihedrals.results`.
        """        

        bins = kwargs.pop('bins', 361)

        results = AttributeDict()

        # get graphs that were produced by g_angle 
        for name, f in self.parameters.filenames.items():
            try:
                results[name] = XVG(f)
            except IOError:
                pass    # either not computed (yet) or some failure

        # compute individual distributions
        ts = results['timeseries'].array    # ts[0] = time, ts[1] = avg
        dih = ts[2:]

        phi_range = (-180., 180.)

        Ndih = len(dih)
        p = Ndih * [None]  # histograms (prob. distributions), one for each dihedral i
        for i in xrange(Ndih):
            phis = dih[i]
            p[i],e = numpy.histogram(phis, bins=bins, range=phi_range, normed=True, new=True)

        P = numpy.array(p)
        phi = 0.5*(e[:-1]+e[1:])   # midpoints of bin edges
        distributions = numpy.concatenate((phi[numpy.newaxis, :], P))  # phi, P[0], P[1], ...

        xvg = XVG()
        xvg.set(distributions)
        xvg.write(self.parameters.filenames['distributions'])
        results['distributions'] = xvg
        del xvg
        
        # compute PMF (from individual distributions)
        W = -numpy.log(P)                      # W(phi)/kT = -ln P
        W -= W.min(axis=1)[:, numpy.newaxis]   # minimum at 0 kT
        pmf = numpy.concatenate((phi[numpy.newaxis, :], W), axis=0)
        xvg = XVG()
        xvg.set(pmf)
        xvg.write(self.parameters.filenames['PMF'])
        results['PMF'] = xvg
        
        self.results = results
        return results

    def plot(self, **kwargs):
        """Plot all results in one graph, labelled by the result keys.

        :Keywords:
           figure
               - ``True``: save figures in the given formats
               - "name.ext": save figure under this filename (``ext`` -> format)
               - ``False``: only show on screen
           formats : sequence
               sequence of all formats that should be saved [('png', 'pdf')]
           with_legend : bool
               add legend from *labels* to graphs
           plotargs    
               keyword arguments for pylab.plot()
        """

        from pylab import subplot, xlabel, ylabel, legend
        figure = kwargs.pop('figure', False)
        extensions = kwargs.pop('formats', ('pdf','png'))

        kwargs.setdefault('linestyle', '-')
        kwargs.setdefault('linewidth', 3)
        with_legend = kwargs.pop('with_legend', True)

        # data
        try:
            P = self.results['distributions']
            W = self.results['PMF']
        except KeyError:
            raise MissingDataError("No post-processed data available: "
                                   "Run the analysis() method first.")

        # plot probability and PMF
        subplot(211)
        P.plot(**kwargs)
        subplot(212)
        W.plot(**kwargs)
        xlabel('dihedral angle $\phi/\degree$')
        ylabel(r'PMF  $\mathcal{W}/kT$')

        # use legend labels from init
        if with_legend and len(self.parameters.labels) > 0:
            if 'columns' in kwargs:
                # sneaky option for XVG.plot: selects columns to plot, eg [0,1, 3] so we
                # need to get the correct labels, namely [0, 2]
                c = numpy.asarray(kwargs['columns'])
                not_zero = (c != 0)
                c = c[not_zero]  # filter 0th column (phi) which does not have a label
                c -= 1         # it's now an index into labels
                if numpy.all(not_zero):
                    # special case when 0th column was not included (eg correlation plots):
                    # we have to discard the first label (and *should* change the xlabel...)
                    xlabel(self.parameters.labels[c[0]])
                    c = c[1:]
                _labels = tuple([self.parameters.labels[i] for i in c])
            else:
                _labels = tuple(self.parameters.labels)
            subplot(211)
            legend(_labels, loc='best')

        # write graphics file(s)
        if figure is True:
            for ext in extensions:
                self.savefig(ext=ext)
        elif figure:
            self.savefig(filename=figure)



# Public classes that register the worker classes
#------------------------------------------------

class Dihedrals(Plugin):
    """*Dihedrals* plugin.
    
    .. class:: Dihedrals(dihedrals[, labels[,name[, simulation]]])
    
    :Keywords:
        *dihedrals*        
            list of tuples; each tuple contains atom indices that define the dihedral.            
        *labels*
            optional list of labels for the dihedrals. Must have as many entries as
            *dihedrals*.
        *name* : string
            plugin name (used to access it)
        *simulation* : instance
            The :class:`gromacs.analysis.Simulation` instance that owns the plugin.

    """
    worker_class = _Dihedrals



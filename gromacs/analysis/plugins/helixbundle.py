# plugin: helixbundle.py
# Copyright (c) 2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
Helix bundle analysis
=====================

Analysis of helix bundles with :func:`gromacs.g_bundle`.

.. SeeAlso:: HELANAL does some things better than g_bundle, in
              particular it can find kinks automatically.

Plugin class
------------

.. autoclass:: HelixBundle
   :members: worker_class
   :undoc-members:

Example
~~~~~~~

.. Warning:: Note that the helices should be sorted in ascending order
             because all indices are sorted that way. If they are in
             mixed order then the correspondence between the name
             column and the helices is lost. This is a shortcoming of
             :program:`make_ndx` that is not easy to rectify.

Define helices in a reST table (use column names as shown!), set up
the plugin, and add it to a simulation instance::

   helixtable = '''
        Table[helices]: helices top = N-term, bot = C-term; kink
        =====  ================ ================ ===============
        name   top              bottom           kink
        =====  ================ ================ ===============
        TM1    G28  P29  F30    T53  S54  S55    I41  Q42  V43
        TM2    Q57  V58  W59    I84  R85  W86    F78  T79  Q80
        TM3    R101 G102 S103   L135 L136 T137   F116 W117 F118
        TM4    L142 P143 L144   T157 F158 Y159   F149 G150 A151
        TM6    F209 S210 T211   D229 I230 V231   G219 W220 I221
        TM7    E243 G244 Q245   L277 V278 G279   V261 P262 A263
        TM8    P297 M298 A299   C327 S328 T329   N314 P315 A316
        TM9    F336 K337 T338   L348 L349 M350   V342 S343 A344
        =====  ================ ================ ===============
        '''
   hb = HelixBundle(helixtable=helixtable, offset=-9)
   S = Simulation(....)
   S.add_plugin(hb)


Worker class
------------

The worker class performs the analysis.

.. autoclass:: _HelixBundle
   :members:


"""
from __future__ import with_statement

__docformat__ = "restructuredtext en"

import os.path
import warnings

import recsql

import gromacs
from gromacs.utilities import AttributeDict
from gromacs.analysis.core import Worker, Plugin

import logging
logger = logging.getLogger('gromacs.analysis.plugins.helixbundle')

# Worker classes that are registered via Plugins (see below)
# ----------------------------------------------------------
# These must be defined before the plugins.

class _HelixBundle(Worker):
    """HelixBundle worker class."""

    def __init__(self,**kwargs):
        """Set up HelixBundle analysis.

        :Keywords:
           *helixtable*
               reST table with columns "name", "top", "bottom",
               "kink"; see :mod:`gromacs.analysis.plugins.helixbundle`
               for details
           *offset*
               add the *offset* to the residue numbers in *helixtable* [0]
           *helixndx*
               provide a index file with the appropriate groups
               instead of the table; also requires *na*
           *na*
               number of helices
           *with_kinks*
               take kinks into account [True]
           *name*
               plugin name [HelixBundle]
           *simulation*
               The :class:`gromacs.analysis.Simulation` instance that
               owns the plugin [None]
        """
        helixndx = kwargs.pop('helixndx', None)
        helixtable = kwargs.pop('helixtable', None)
        na = kwargs.pop('na', None)
        offset = kwargs.pop('offset', 0)
        with_kinks = kwargs.pop('with_kinks', True)

        if helixndx is None:
            if helixtable is None:
                raise ValueError("HelixBundle: requires a table *helixtable* with residues "
                                 "that indicate the helix boundaries and kink or a index "
                                 "file *helixndx* appropriate for g_bundle. See the "
                                 "docs for details.")
            self.helices = recsql.rest_table.Table2array(helixtable, autoconvert=True).recarray()
            na = len(self.helices)
        elif na is None:
            raise ValueError("When using *helixndx* one must also provide the total number of "
                             "helices *na*. See g_bundle docs for details.")

        super(_HelixBundle, self).__init__(**kwargs)
        
        self.parameters.na = na
        self.parameters.offset = offset
        self.parameters.with_kinks = with_kinks

        if not self.simulation is None:
            self._register_hook()

    def _register_hook(self, **kwargs):
        """Run when registering; requires simulation."""

        super(_HelixBundle, self)._register_hook(**kwargs)
        assert not self.simulation is None

        self.helixndx = self.plugindir('helices.ndx')     # special index for g_bundle

        self.parameters.filenames = {                     # result xvg files
            'length': self.plugindir('len.xvg'),
            'distance': self.plugindir('dist.xvg'),
            'z': self.plugindir('z.xvg'),
            'tilt': self.plugindir('tilt.xvg'),
            'tilt_radial': self.plugindir('tilt_radial.xvg'),
            'tilt_lateral': self.plugindir('tilt_lateral.xvg'),
            }
        if self.parameters.with_kinks:
            self.parameters.filenames.update(
                {'kink': self.plugindir('kink.xvg'),
                 'kink_radial': self.plugindir('kink_radial.xvg'),
                 'kink_lateral': self.plugindir('kink_lateral.xvg'),
                 })

        # default filename for the plots -- not used
        self.parameters.fignames = {
            'z': self.figdir('z'),
            'tilt': self.figdir('tilt'),
            'kink': self.figdir('kink'),
            }
            
    def make_index(self, force=None):
        """Build g_bundle index file from a record array.

        :Keywords:
          *force*
              - ``None`` does nothing if the index file already exists
              - ``False`` raises exception if it exists
              - ``True`` always rebuild index file

        Format of the index table:
        - must contain columns name, top, bottom, kink
        - a row corresponds to one helis
        - name contains the name of the helix e.g. 'TM5'
        - a entry is a *string* of residues which is split on white space; 
          entries in the same column must have the same number of residues 
          (limitation of g_bundle)

        """

        if self.check_file_exists(self.helixndx, resolve="indicate", force=force):
            logger.warning("make_index(): The index file %r already exists; "
                           "use force = True to rebuild", self.helixndx)
            return self.helixndx

        logger.info("Making index file for %d helices %r", len(self.helices), self.helixndx)
        logger.debug("make_index: tpr = %r", self.simulation.tpr)
        logger.debug("make_index: offset = %r", self.parameters.offset)

        # quick'n'crappy...
        tops = []
        bottoms = []
        kinks = []

        for h in self.helices:
            # split entries on white space (Table2array() does not handle it properly yet)
            tops.extend(h.top.split())
            bottoms.extend(h.bottom.split())
            kinks.extend(h.kink.split())

        tpr = self.simulation.tpr
        offset = self.parameters.offset

        # TODO: use tmp files for indices

        Itop = gromacs.cbook.IndexBuilder(tpr, tops, offset=offset, out_ndx='top.ndx')
        name, topndx = Itop.combine(name_all='top', defaultgroups=False)

        Ibot = gromacs.cbook.IndexBuilder(tpr, bottoms, offset=offset, out_ndx='bottom.ndx')
        name, botndx = Ibot.combine(name_all='bottom', defaultgroups=False)

        Ikink = gromacs.cbook.IndexBuilder(tpr, kinks, offset=offset, out_ndx='kinks.ndx')
        name, kinkndx = Ikink.combine(name_all='kink', defaultgroups=False)

        Iall = gromacs.cbook.IndexBuilder(tpr, ndx=[topndx, botndx, kinkndx])
        Iall.cat(self.helixndx)

        return self.helixndx



    def run(self, force=None, **gmxargs):
        """Analyze trajectory and write HelixBundle files.

        :Arguments:
          - *force*: ``True`` does analysis and overwrites existing files
          - *gmxargs*: additional keyword arguments for :func:`gromacs.g_bundle` 

        .. Note:: The plugin default is *z* = ``True``, i.e. the tilt is computed
                  relative to the box z-axis.
        """
        gmxargs.setdefault('z', True)
        gmxargs['na'] = self.parameters.na

        if self.check_file_exists(self.parameters.filenames['tilt'], resolve='warning', force=force):
            return

        self.make_index(force=force)

        logger.info("Analyzing HelixBundle (with_kinks=%r)...", self.parameters.with_kinks)
        f = self.parameters.filenames
        if self.parameters.with_kinks:
            gromacs.g_bundle(s=self.simulation.tpr, f=self.simulation.xtc, n=self.helixndx,
                             ol=f['length'], od=f['distance'], oz=f['z'], 
                             ot=f['tilt'], otr=f['tilt_radial'], otl=f['tilt_lateral'],
                             ok=f['kink'], okr=f['kink_radial'], okl=f['kink_lateral'],
                             input=['top','bottom', 'kink'],
                             **gmxargs)
        else:
            gromacs.g_bundle(s=self.simulation.tpr, f=self.simulation.xtc, n=self.helixndx,
                             ol=f['length'], od=f['distance'], oz=f['z'], 
                             ot=f['tilt'], otr=f['tilt_radial'], otl=f['tilt_lateral'],
                             input=['top','bottom'],
                             **gmxargs)


    def analyze(self,**kwargs):
        """Collect output xvg files as :class:`gromacs.formats.XVG` objects.

        :Returns:  a dictionary of the results and also sets ``self.results``.
        """        
        from gromacs.formats import XVG

        logger.info("Preparing HelixBundle graphs as XVG objects.")        
        results = AttributeDict( (k, XVG(fn)) for k,fn in self.parameters.filenames.items() )
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
           plotargs    
               keyword arguments for pylab.plot()
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
        if figure is True:
            for ext in extensions:
                self.savefig(ext=ext)
        elif figure:
            self.savefig(filename=figure)

    


# Public classes that register the worker classes
#------------------------------------------------

class HelixBundle(Plugin):
    """*HelixBundle* plugin.

    :func:`gromacs.g_bundle` helix analysis

    .. class:: HelixBundle([helixtable, offset, with_kinks, [name[, simulation]]])
    
    :Arguments:
       *helixtable*
           reST table with columns "name", "top", "bottom",
           "kink"; see :mod:`gromacs.analysis.plugins.helixbundle`
           for details
       *offset*
           add the *offset* to the residue numbers in *helixtable* [0]
       *helixndx*
           provide a index file with the appropriate groups
           instead of the table; also requires *na*
       *na*
           number of helices
       *with_kinks*
           take kinks into account [True]
       *name*
           plugin name [HelixBundle]
       *simulation*
           The :class:`gromacs.analysis.Simulation` instance that owns
           the plugin. [None]
    """
    worker_class = _HelixBundle



# $Id$
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

"""
Various terms from the energy (edr) file
========================================

Simplified invocation of :func:`gromacs.g_energy`; only a few selected
terms are plotted.


Plugin class
------------

.. autoclass:: Energy
   :members: worker_class
   :undoc-members:

Worker class
------------

The worker class performs the analysis.

.. autoclass:: _Energy
   :members:


"""
from __future__ import with_statement

__docformat__ = "restructuredtext en"

import os.path
import warnings

import gromacs
from gromacs.utilities import AttributeDict
from gromacs.analysis.core import Worker, Plugin

import logging
logger = logging.getLogger('gromacs.analysis.plugins.energy')

# Worker classes that are registered via Plugins (see below)
# ----------------------------------------------------------
# These must be defined before the plugins.

class _Energy(Worker):
    """Energy worker class.

    **Gromacs 4-git energy terms**

    Select the terms you want from the following list by selecting
    either (part of) the name or the number or a combination.  End
    your selection with an empty line or a zero ::

      1  Bond             2  Angle            3  Proper-Dih.      4  Ryckaert-Bell.
      5  Improper-Dih.    6  LJ-14            7  Coulomb-14       8  LJ-(SR)       
      9  Disper.-corr.   10  Coulomb-(SR)    11  Coul.-recip.    12  Potential     
     13  Kinetic-En.     14  Total-Energy    15  Temperature     16  Pressure-(bar)
     17  Cons.-rmsd-()   18  Box-X           19  Box-Y           20  Box-Z         
     21  Volume          22  Density-(SI)    23  pV              24  Vir-XX        
     25  Vir-XY          26  Vir-XZ          27  Vir-YX          28  Vir-YY        
     29  Vir-YZ          30  Vir-ZX          31  Vir-ZY          32  Vir-ZZ        
     33  Pres-XX-(bar)   34  Pres-XY-(bar)   35  Pres-XZ-(bar)   36  Pres-YX-(bar) 
     37  Pres-YY-(bar)   38  Pres-YZ-(bar)   39  Pres-ZX-(bar)   40  Pres-ZY-(bar) 
     41  Pres-ZZ-(bar)   42  #Surf*SurfTen   43  Mu-X            44  Mu-Y          
     45  Mu-Z                                46  Coul-SR:SOLVENT-SOLVENT           
     47  LJ-SR:SOLVENT-SOLVENT               48  Coul-14:SOLVENT-SOLVENT           
     49  LJ-14:SOLVENT-SOLVENT               50  Coul-SR:SOLVENT-LIPIDS            
     51  LJ-SR:SOLVENT-LIPIDS                52  Coul-14:SOLVENT-LIPIDS            
     53  LJ-14:SOLVENT-LIPIDS                54  Coul-SR:SOLVENT-Protein           
     55  LJ-SR:SOLVENT-Protein               56  Coul-14:SOLVENT-Protein           
     57  LJ-14:SOLVENT-Protein               58  Coul-SR:LIPIDS-LIPIDS             
     59  LJ-SR:LIPIDS-LIPIDS                 60  Coul-14:LIPIDS-LIPIDS             
     61  LJ-14:LIPIDS-LIPIDS                 62  Coul-SR:LIPIDS-Protein            
     63  LJ-SR:LIPIDS-Protein                64  Coul-14:LIPIDS-Protein            
     65  LJ-14:LIPIDS-Protein                66  Coul-SR:Protein-Protein           
     67  LJ-SR:Protein-Protein               68  Coul-14:Protein-Protein           
     69  LJ-14:Protein-Protein               70  T-SOLVENT                         
     71  T-LIPIDS        72  T-Protein       73  Lamb-SOLVENT    74  Lamb-LIPIDS   
     75  Lamb-Protein  
    """
    
    terms = ['Potential', 'Kinetic-En.', 'Total-Energy',
             'Temperature', 'Pressure',
             'Volume', 'Box-X', 'Box-Y', 'Box-Z',
             ]

    def __init__(self,**kwargs):
        """Set up Energy analysis.

        This is the worker class; this is where all the real analysis is done.
        """
        super(_Energy, self).__init__(**kwargs)
        if not self.simulation is None:
            self._register_hook()

    def _register_hook(self, **kwargs):
        """Run when registering; requires simulation."""

        super(_Energy, self)._register_hook(**kwargs)
        assert not self.simulation is None

        self.parameters.filenames = {
            'Energy': self.plugindir('energy.xvg'),
            }
        # default filename for the plot
        self.parameters.figname = self.figdir('energy')

    def run(self, force=False, **gmxargs):
        """Analyze trajectory and write energy files.

        :Arguments:
          - *force*: do analysis and overwrite existing files
          - *gmxargs*: additional keyword arguments for :func:`gromacs.g_energy`
        """
        if not self.check_file_exists(self.parameters.filenames['Energy'], resolve='warning') or force:
            logger.info("Analyzing energy file...")
            gromacs.g_energy(s=self.simulation.tpr, f=self.simulation.edr,
                             o=self.parameters.filenames['Energy'],
                             input=self.terms + [0], **gmxargs)

    def analyze(self,**kwargs):
        """Collect output xvg files as :class:`gromacs.formats.XVG` objects.

        :Returns:  a dictionary of the results and also sets ``self.results``.
        """        
        from gromacs.formats import XVG

        logger.info("Preparing Energy graphs as XVG objects.")
        results = AttributeDict(Energy=XVG(self.parameters.filenames['Energy']))
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

class Energy(Plugin):
    """*Energy* plugin.
    
    Analysis of terms in the Gromacs energy (edr) file.

    .. class:: Energy([name[, simulation]])

     The following energy terms are extracted from the edr file: 

       * 'Potential'
       * 'Kinetic-En.'
       * 'Total-Energy'
       * 'Temperature'
       * 'Pressure'
       * 'Volume'
       * 'Box-X', 'Box-Y', 'Box-Z'

    The list of terms corresponds to input values for :func:`gromacs.g_energy`
    and is defined in the attribute :attr:`_Energy.terms`.
    
    :Arguments:
        *name* : string
            plugin name (used to access it)
        *simulation* : instance
            The :class:`gromacs.analysis.Simulation` instance that owns the plugin.

    """
    worker_class = _Energy



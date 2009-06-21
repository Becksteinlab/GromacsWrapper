# $Id$
"""Setting up a Gromacs MD run.

Individual steps such as solvating a structure or energy minimization
are set up in individual directories.

"""

from __future__ import with_statement
import os
import errno
import re
import shutil
from contextlib import contextmanager

from pkg_resources import resource_filename


import gromacs
import gromacs.cbook

@contextmanager
def in_dir(directory):
    """Context manager to execute a code block in a directory.

    * directory is created if it does not exist
    * at the end or after and exception code always returns to
      the directory that was the current directory before entering
      the block
    """
    startdir = os.getcwd()
    try:
        try:
            os.chdir(directory)
        except OSError, err:
            if err.errno == errno.ENOENT:
                os.makedirs(directory)
                os.chdir(directory)
            else:
                raise
        yield os.getcwd()
    finally:
        os.chdir(startdir)

# TODO:
# - should be part of a class so that we can store the topology etc
# - full logging would be nice (for provenance)

def topology(struct=None, protein='protein',
             ff='oplsaa', watermodel='TIP4P', 
             top='system.top',  dirname='top'):
    """Build Gromacs topology files from pdb.
    
    INCOMPLETE
    """
    raise NotImplemented('sorry, this needs to be done manually at the moment')

    structure = os.path.realpath(struct)
    topology = os.path.realpath(os.path.join(dirname, top))
    
    new_struct = protein + '.pdb'
    posres = protein + '_posres.itp'

    with in_dir(dirname):
        gromacs.pdb2gmx(ff=ff, water=watermodel, f=structure,
                        o=new_struct, p='tmp', i=posres)
        # need some editing  AKeco_tmp --> system.top and AKeco.itp

    return topology, new_struct

def solvate(struct='top/protein.pdb', top='top/system.top',
            distance=0.9, boxtype='dodecahedron',
            concentration=0, cation='NA+', anion='CL-',
            dirname='solvate'):
    """Put protein into box, add water, add ions."""

    structure = os.path.realpath(struct)
    topology = os.path.realpath(top)

    # should scrub topology but hard to know what to scrub; for
    # instance, some ions or waters could be part of the crystal structure
    
    with in_dir(dirname):
        gromacs.editconf(f=structure, o='boxed.gro', bt=boxtype, d=distance)
        gromacs.genbox(p=topology, cp='boxed.gro', cs='tip4p', o='solvated.gro')

        with open('none.mdp','w') as mdp:
            print >> mdp, '; empty'

        qtot = gromacs.cbook.grompp_qtot(f='none.mdp', o='topol.tpr', c='solvated.gro',
                                         p=topology, stdout=False)
        print "After solvation: total charge qtot = %(qtot)r" % vars()        

        # TODO:
        # - get number of waters <---- ! hold up (grep topol...)
        # - calculate numbers from concentration
        # target concentration of free ions: c = 100mM  (c_water = 55M)
        # N = N_water * c/c_water
        # add ions for concentration to the counter ions (counter ions are less free)
        if concentration != 0:
            raise NotImplementedError('only conc = 0 working at moment')

        # neutralize (or try -neutral switch of genion)
        n_cation = n_anion = 0
        if qtot > 0:
            n_anion = int(abs(qtot))
        elif qtot < 0:
            n_cation = int(abs(qtot))
        gromacs.genion(s='topol.tpr', o='ionized.gro', p=topology,
                       pname=cation, nname=anion, np=n_cation, nn=n_anion,
                       input='SOL')

        qtot = gromacs.cbook.grompp_qtot(f='none.mdp', o='ionized.tpr', c='ionized.gro',
                                         p=topology, stdout=False)
        gromacs.cbook.trj_compact(f='ionized.gro', s='ionized.tpr', o='compact.pdb')
        return qtot

def energy_minimize(mdp=resource_filename('templates/em.mdp'),
                    struct='solvate/ionized.gro', top='top/system.top', dirname='em'):
    """Energy minimize the system."""
    structure = os.path.realpath(struct)
    topology = os.path.realpath(top)
    mdp = os.path.realpath(mdp)

    tpr = 'em.tpr'
    with in_dir(dirname):
        gromacs.grompp(f=mdp, o=tpr, c=structure, p=topology, maxwarn=1)
        # TODO: not clear yet how to run as MPI with mpiexec & friends
        # TODO: fall back to mdrun wif no double precision binary
        gromacs.mdrun_d('v', stepout=10, deffnm='em', c='em.pdb')   # or em.pdb ? (box??)
        # em.gro --> gives 'Bad box in file em.gro' warning ??


def _setup_MD(dirname,
              deffnm='md', mdp=resource_filename('templates/md.mdp'), struct=None,
              top='top/system.top', ndx=None,
              sge=resource_filename('templates/deathspud.sge'),
              dt=0.002, runtime=1e3, **mdp_kwargs):
    structure = os.path.realpath(struct)
    topology = os.path.realpath(top)
    mdp_template = os.path.realpath(mdp)
    sge_template = os.path.realpath(sge)

    nsteps = int(float(runtime)/float(dt))

    mdp = deffnm + '.mdp'
    tpr = deffnm + '.tpr'

    mdp_parameters = {'nsteps':nsteps, 'dt':dt}
    mdp_parameters.update(mdp_kwargs)
    
    with in_dir(dirname):
        # do I need an index file?
        gromacs.cbook.edit_mdp(mdp_template, new_mdp=mdp, 
                               **mdp_parameters)
        gromacs.grompp(f=mdp, p=topology, c=structure, n=ndx, o=tpr)
        # edit sge script manually
        # TODO: edit DEFFNM in sge template
        shutil.copy(sge_template, os.path.curdir)        

    print "All files set up for a run time of %(runtime)g ps "\
        "(dt=%(dt)g, nsteps=%(nsteps)g)" % vars()
    return {'dirname':dirname, 'tpr':tpr}


def MD_restrained(dirname='MD_POSRES', **kwargs):
    kwargs.setdefault('struct', 'em/em.pdb')
    kwargs['define'] = '-DPOSRES'
    return _setup_MD(dirname, **kwargs)

def MD(dirname='MD', **kwargs):
    kwargs.setdefault('struct', 'md_posres/md_posres.pdb')
    return _setup_MD(dirname, **kwargs)



# TODO: autorun (qv MI lipids/setup.pl)

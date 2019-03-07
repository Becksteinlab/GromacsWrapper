# GromacsWrapper: setup.py
"""
:mod:`gromacs.setup` -- Setting up a Gromacs MD run
===================================================

Individual steps such as solvating a structure or energy minimization
are set up in individual directories. For energy minimization one
should supply appropriate mdp run input files; otherwise example
templates are used.

.. warning::

   You **must** check all simulation parameters for yourself. Do not rely on
   any defaults provided here. The scripts provided here are provided under the
   assumption that you know what you are doing and you just want to automate
   the boring parts of the process.


User functions
--------------

The individual steps of setting up a simple MD simulation are broken down in a
sequence of functions that depend on the previous step(s):

  :func:`topology`
        generate initial topology file (limited functionality, might require
        manual setup)
  :func:`solvate`
        solvate globular protein and add ions to neutralize
  :func:`energy_minimize`
        set up energy minimization and run it (using ``mdrun_d``)
  :func:`em_schedule`
        set up and run multiple energy minimizations one after another (as an
        alternative to the simple single energy minimization provided by
        :func:`energy_minimize`)
  :func:`MD_restrained`
        set up restrained MD
  :func:`MD`
        set up equilibrium MD

Each function uses its own working directory (set with the ``dirname`` keyword
argument, but it should be safe and convenient to use the defaults). Other
arguments assume the default locations so typically not much should have to be
set manually.

One can supply non-standard itp files in the topology directory. In
some cases one does not use the :func:`topology` function at all but
sets up the topology manually. In this case it is safest to call the
topology directory ``top`` and make sure that it contains all relevant
top, itp, and pdb files.


Example
-------

Run a single protein in a dodecahedral box of SPC water molecules and
use the GROMOS96 G43a1 force field. We start with the structure in
``protein.pdb``::

  from gromacs.setup import *
  f1 = topology(protein='MyProtein', struct='protein.pdb', ff='G43a1', water='spc', force=True, ignh=True)

Each function returns "interesting" new files in a dictionary in such
a away that it can often be used as input for the next function in the
chain (although in most cases one can get away with the defaults of
the keyword arguments)::

  f2 = solvate(**f1)
  f3 = energy_minimize(**f2)

Now prepare input for a MD run with restraints on the protein::

  MD_restrained(**f3)

Use the files in the directory to run the simulation locally or on a
cluster. You can provide your own template for a queuing system
submission script; see the source code for details.

Once the restraint run has completed, use the last frame as input for
the equilibrium MD::

  MD(struct='MD_POSRES/md.gro', runtime=1e5)

Run the resulting tpr file on a cluster.


User functions
--------------

The following functions are provided for the user:

.. autofunction:: topology
.. autofunction:: solvate
.. autofunction:: energy_minimize
.. autofunction:: em_schedule
.. autofunction:: MD_restrained
.. autofunction:: MD

Helper functions
----------------

The following functions are used under the hood and are mainly useful when
writing extensions to the module.

.. autofunction:: make_main_index
.. autofunction:: check_mdpargs
.. autofunction:: get_lipid_vdwradii
.. autofunction:: _setup_MD

Defined constants:

.. autodata:: CONC_WATER
.. autodata:: vdw_lipid_resnames
.. autodata:: vdw_lipid_atom_radii

"""

from __future__ import absolute_import, with_statement

__docformat__ = "restructuredtext en"

import os
import errno
import re
import shutil
import warnings

import logging
logger = logging.getLogger('gromacs.setup')


import gromacs
from . import config
from .exceptions import (GromacsError, GromacsFailureWarning, GromacsValueWarning,
                         AutoCorrectionWarning, BadParameterWarning, UsageWarning,
                         MissingDataError)
from . import run
from . import cbook
from . import qsub
from . import utilities
from .utilities import in_dir, realpath, Timedelta, asiterable, firstof


#: Concentration of water at standard conditions in mol/L.
#: Density at 25 degrees C and 1 atmosphere pressure: rho = 997.0480 g/L.
#: Molecular weight: M = 18.015 g/mol.
#: c = n/V = m/(V*M) = rho/M  = 55.345 mol/L.
CONC_WATER = 55.345

# XXX: This is not used anywhere at the moment:
# parse this table for a usable data structure (and then put it directly in the docs)
recommended_mdp_table = """\
Table: recommended mdp parameters for different FF
==========   =========  ================
mdp          GROMOS     OPLS-AA
==========   =========  ================
rvdw         1.4        1.0
rlist        1.4 ?      1.0
==========   =========  ================
"""


trj_compact_main = gromacs.tools.Trjconv(ur='compact', center=True, boxcenter='tric', pbc='mol',
                                         input=('__main__', 'system'))
# trj_compact_main.__doc__ += "Returns a compact representation of the system centered on the __main__ group"

# TODO:
# - should be part of a class so that we can store the topology etc !!!
#   and also store mainselection

def topology(struct=None, protein='protein',
             top='system.top',  dirname='top',
             posres="posres.itp",
             ff="oplsaa", water="tip4p",
             **pdb2gmx_args):
    """Build Gromacs topology files from pdb.

    :Keywords:
       *struct*
           input structure (**required**)
       *protein*
           name of the output files
       *top*
           name of the topology file
       *dirname*
           directory in which the new topology will be stored
       *ff*
           force field (string understood by ``pdb2gmx``); default
           "oplsaa"
       *water*
           water model (string), default "tip4p"
       *pdb2gmxargs*
           other arguments for ``pdb2gmx``

    .. note::
       At the moment this function simply runs ``pdb2gmx`` and uses
       the resulting topology file directly. If you want to create
       more complicated topologies and maybe also use additional itp
       files or make a protein itp file then you will have to do this
       manually.
    """

    structure = realpath(struct)

    new_struct = protein + '.pdb'
    if posres is None:
        posres = protein + '_posres.itp'

    pdb2gmx_args.update({'f': structure, 'o': new_struct, 'p': top, 'i': posres,
                         'ff': ff, 'water': water})

    with in_dir(dirname):
        logger.info("[{dirname!s}] Building topology {top!r} from struct = {struct!r}".format(**vars()))
        # perhaps parse output from pdb2gmx 4.5.x to get the names of the chain itp files?
        gromacs.pdb2gmx(**pdb2gmx_args)
    return { \
            'top': realpath(dirname, top), \
            'struct': realpath(dirname, new_struct), \
            'posres' : realpath(dirname, posres) }

def make_main_index(struct, selection='"Protein"', ndx='main.ndx', oldndx=None):
    """Make index file with the special groups.

    This routine adds the group __main__ and the group __environment__
    to the end of the index file. __main__ contains what the user
    defines as the *central* and *most important* parts of the
    system. __environment__ is everything else.

    The template mdp file, for instance, uses these two groups for T-coupling.

    These groups are mainly useful if the default groups "Protein" and "Non-Protein"
    are not appropriate. By using symbolic names such as __main__ one
    can keep scripts more general.

    :Returns:
      *groups* is a list of dictionaries that describe the index groups. See
      :func:`gromacs.cbook.parse_ndxlist` for details.

    :Arguments:
      *struct* : filename
        structure (tpr, pdb, gro)
      *selection* : string
        is a ``make_ndx`` command such as ``"Protein"`` or ``r DRG`` which
        determines what is considered the main group for centering etc. It is
        passed directly to ``make_ndx``.
      *ndx* : string
         name of the final index file
      *oldndx* : string
         name of index file that should be used as a basis; if None
         then the ``make_ndx`` default groups are used.

    This routine is very dumb at the moment; maybe some heuristics will be
    added later as could be other symbolic groups such as __membrane__.
    """

    logger.info("Building the main index file {ndx!r}...".format(**vars()))

    # pass 1: select
    # get a list of groups
    # need the first "" to get make_ndx to spit out the group list.
    _,out,_ = gromacs.make_ndx(f=struct, n=oldndx, o=ndx, stdout=False,
                                      input=("", "q"))
    groups = cbook.parse_ndxlist(out)

    # find the matching groups,
    # there is a nasty bug in GROMACS where make_ndx may have multiple
    # groups, which caused the previous approach to fail big time.
    # this is a work around the make_ndx bug.
    # striping the "" allows compatibility with existing make_ndx selection commands.
    selection = selection.strip("\"")

    selected_groups = [g for g in groups if g['name'].lower() == selection.lower()]

    if len(selected_groups) > 1:
        logging.warn("make_ndx created duplicated groups, performing work around")

    if len(selected_groups) <= 0:
        msg = "no groups found for selection {0}, available groups are {1}".format(selection, groups)
        logging.error(msg)
        raise ValueError(msg)

    # Found at least one matching group, we're OK

    # index of last group
    last = len(groups) - 1
    assert last == groups[-1]['nr']

    group = selected_groups[0]

    # pass 2:
    # 1) last group is __main__
    # 2) __environment__ is everything else (eg SOL, ions, ...)
    _,out,_ = gromacs.make_ndx(f=struct, n=ndx, o=ndx,
                                      stdout=False,
                                             # make copy selected group, this now has index last + 1
                                      input=("{0}".format(group['nr']),
                                             # rename this to __main__
                                             "name {0} __main__".format(last+1),
                                             # make a complement to this group, it get index last + 2
                                             "! \"__main__\"",
                                             # rename this to __environment__
                                             "name {0} __environment__".format(last+2),
                                             # list the groups
                                             "",
                                             # quit
                                             "q"))
    return cbook.parse_ndxlist(out)


#: Hard-coded lipid residue names for a ``vdwradii.dat`` file. Use together with
#: :data:`~gromacs.setup.vdw_lipid_atom_radii` in :func:`~gromacs.setup.get_lipid_vdwradii`.
vdw_lipid_resnames = ["POPC", "POPE", "POPG", "DOPC", "DPPC", "DLPC", "DMPC", "DPPG"]
#: Increased atom radii for lipid atoms; these are simply the standard values from
#: ``GMXLIB/vdwradii.dat`` increased by 0.1 nm (C) or 0.05 nm (N, O, H).
vdw_lipid_atom_radii = {'C': 0.25, 'N': 0.16, 'O': 0.155, 'H': 0.09}

def get_lipid_vdwradii(outdir=os.path.curdir, libdir=None):
    """Find vdwradii.dat and add special entries for lipids.

    See :data:`gromacs.setup.vdw_lipid_resnames` for lipid
    resnames. Add more if necessary.
    """
    vdwradii_dat = os.path.join(outdir, "vdwradii.dat")

    if libdir is not None:
        filename = os.path.join(libdir, 'vdwradii.dat')  # canonical name
        if not os.path.exists(filename):
            msg = 'No VDW database file found in {filename!r}.'.format(**vars())
            logger.exception(msg)
            raise OSError(msg, errno.ENOENT)
    else:
        try:
            filename = os.path.join(os.environ['GMXLIB'], 'vdwradii.dat')
        except KeyError:
            try:
                filename = os.path.join(os.environ['GMXDATA'], 'top', 'vdwradii.dat')
            except KeyError:
                msg = "Cannot find vdwradii.dat. Set GMXLIB (point to 'top') or GMXDATA ('share/gromacs')."
                logger.exception(msg)
                raise OSError(msg, errno.ENOENT)
    if not os.path.exists(filename):
        msg = "Cannot find {filename!r}; something is wrong with the Gromacs installation.".format(**vars())
        logger.exception(msg, errno.ENOENT)
        raise OSError(msg)

    # make sure to catch 3 and 4 letter resnames
    patterns = vdw_lipid_resnames + list({x[:3] for x in vdw_lipid_resnames})
    # TODO: should do a tempfile...
    with open(vdwradii_dat, 'w') as outfile:
        # write lipid stuff before general
        outfile.write('; Special larger vdw radii for solvating lipid membranes\n')
        for resname in patterns:
            for atom,radius in vdw_lipid_atom_radii.items():
                outfile.write('{resname:4!s} {atom:<5!s} {radius:5.3f}\n'.format(**vars()))
        with open(filename, 'r') as infile:
            for line in infile:
                outfile.write(line)
    logger.debug('Created lipid vdW radii file {vdwradii_dat!r}.'.format(**vars()))
    return realpath(vdwradii_dat)

def solvate_sol(struct='top/protein.pdb', top='top/system.top',
                distance=0.9, boxtype='dodecahedron',
                water='tip4p', solvent_name='SOL', with_membrane=False,
                dirname='solvate',
                **kwargs):
    structure = realpath(struct)
    topology = realpath(top)

    # arguments for editconf that we honour
    editconf_keywords = ["box", "bt", "angles", "c", "center", "aligncenter",
                         "align", "translate", "rotate", "princ"]
    editconf_kwargs = dict((k,kwargs.pop(k,None)) for k in editconf_keywords)
    editconf_boxtypes = ["triclinic", "cubic", "dodecahedron", "octahedron", None]

    # needed for topology scrubbing
    scrubber_kwargs = {'marker': kwargs.pop('marker',None)}

    # sanity checks and argument dependencies
    bt = editconf_kwargs.pop('bt')
    boxtype = bt if bt else boxtype   # bt takes precedence over boxtype
    if not boxtype in editconf_boxtypes:
        msg = "Unsupported boxtype {boxtype!r}: Only {boxtypes!r} are possible.".format(**vars())
        logger.error(msg)
        raise ValueError(msg)
    if editconf_kwargs['box']:
        distance = None    # if box is set then user knows what she is doing...

    if water.lower() in ('spc', 'spce'):
        water = 'spc216'
    elif water.lower() == 'tip3p':
        water = 'spc216'
        logger.warning("TIP3P water model selected: using SPC equilibrated box "
                       "for initial solvation because it is a reasonable starting point "
                       "for any 3-point model. EQUILIBRATE THOROUGHLY!")

    # clean topology (if user added the marker; the default marker is
    # ; Gromacs auto-generated entries follow:
    n_removed = cbook.remove_molecules_from_topology(topology, **scrubber_kwargs)

    with in_dir(dirname):
        logger.info("[{dirname!s}] Solvating with water {water!r}...".format(**vars()))
        if boxtype is None:
            hasBox = False
            ext = os.path.splitext(structure)[1]
            if ext == '.gro':
                hasBox = True
            elif ext == '.pdb':
                with open(structure) as struct:
                    for line in struct:
                        if line.startswith('CRYST'):
                            hasBox = True
                            break
            if not hasBox:
                msg = "No box data in the input structure {structure!r} and boxtype is set to None".format(**vars())
                logger.exception(msg)
                raise MissingDataError(msg)
            distance = boxtype = None   # ensures that editconf just converts
        editconf_kwargs.update({'f': structure, 'o': 'boxed.gro',
                                'bt': boxtype, 'd': distance})
        gromacs.editconf(**editconf_kwargs)

        if with_membrane:
            vdwradii_dat = get_lipid_vdwradii()  # need to clean up afterwards
            logger.info("Using special vdW radii for lipids {0!r}".format(vdw_lipid_resnames))

        try:
            gromacs.genbox(p=topology, cp='boxed.gro', cs=water, o='solvated.gro')
        except:
            if with_membrane:
                # remove so that it's not picked up accidentally
                utilities.unlink_f(vdwradii_dat)
            raise
        logger.info("Solvated system with %s", water)
    return {'struct': realpath(dirname, 'solvated.gro'),}

def solvate_ion(struct='solvated.gro', top='top/system.top',
                concentration=0, cation='NA', anion='CL',
                solvent_name='SOL', ndx='main.ndx',
                mainselection='"Protein"', dirname='solvate',
                **kwargs):
    structure = realpath(struct)
    topology = realpath(top)
    # By default, grompp should not choke on a few warnings because at
    # this stage the user cannot do much about it (can be set to any
    # value but is kept undocumented...)
    grompp_maxwarn = kwargs.pop('maxwarn',10)
    
    # handle additional include directories (kwargs are also modified!)
    mdp_kwargs = cbook.add_mdp_includes(topology, kwargs)

    with in_dir(dirname):
        with open('none.mdp','w') as mdp:
            mdp.write('; empty mdp file\ninclude = {include!s}\nrcoulomb = 1\nrvdw = 1\nrlist = 1\n'.format(**mdp_kwargs))
        qtotgmx = cbook.grompp_qtot(f='none.mdp', o='topol.tpr', c=structure,
                                    p=topology, stdout=False, maxwarn=grompp_maxwarn)
        qtot = round(qtotgmx)
        logger.info("[{dirname!s}] After solvation: total charge qtot = {qtotgmx!r} = {qtot!r}".format(**vars()))

        if concentration != 0:
            logger.info("[{dirname!s}] Adding ions for c = {concentration:f} M...".format(**vars()))
            # target concentration of free ions c ==>
            #    N = N_water * c/c_water
            # add ions for concentration to the counter ions (counter ions are less free)
            #
            # get number of waters (count OW ... works for SPC*, TIP*P water models)
            rc,output,junk = gromacs.make_ndx(f='topol.tpr', o='ow.ndx',
                                              input=('keep 0', 'del 0', 'a OW*', 'name 0 OW', '', 'q'),
                                              stdout=False)
            groups = cbook.parse_ndxlist(output)
            gdict = {g['name']: g for g in groups}   # overkill...
            N_water = gdict['OW']['natoms']                  # ... but dict lookup is nice
            N_ions = int(N_water * concentration/CONC_WATER) # number of monovalents
        else:
            N_ions = 0

        # neutralize (or try -neutral switch of genion???)
        n_cation = n_anion = 0
        if qtot > 0:
            n_anion = int(abs(qtot))
        elif qtot < 0:
            n_cation = int(abs(qtot))

        n_cation += N_ions
        n_anion  += N_ions

        if n_cation != 0 or n_anion != 0:
            # sanity check:
            assert qtot + n_cation - n_anion < 1e-6
            logger.info("[{dirname!s}] Adding n_cation = {n_cation:d} and n_anion = {n_anion:d} ions...".format(**vars()))
            gromacs.genion(s='topol.tpr', o='ionized.gro', p=topology,
                           pname=cation, nname=anion, np=n_cation, nn=n_anion,
                           input=solvent_name)
        else:
            # fake ionized file ... makes it easier to continue without too much fuzz
            try:
                os.unlink('ionized.gro')
            except OSError as err:
                if err.errno != errno.ENOENT:
                    raise
            os.symlink('solvated.gro', 'ionized.gro')

        qtot = cbook.grompp_qtot(f='none.mdp', o='ionized.tpr', c='ionized.gro',
                                 p=topology, stdout=False, maxwarn=grompp_maxwarn)

        if abs(qtot) > 1e-4:
            wmsg = "System has non-zero total charge qtot = {qtot:g} e.".format(**vars())
            warnings.warn(wmsg, category=BadParameterWarning)
            logger.warn(wmsg)

        # make main index
        try:
            make_main_index('ionized.tpr', selection=mainselection, ndx=ndx)
        except GromacsError as err:
            # or should I rather fail here?
            wmsg = "Failed to make main index file %r ... maybe set mainselection='...'.\n"\
                   "The error message was:\n%s\n" % (ndx, str(err))
            logger.warn(wmsg)
            warnings.warn(wmsg, category=GromacsFailureWarning)
        try:
            trj_compact_main(f='ionized.gro', s='ionized.tpr', o='compact.pdb', n=ndx)
        except GromacsError as err:
            wmsg = "Failed to make compact pdb for visualization... pressing on regardless. "\
                   "The error message was:\n%s\n" % str(err)
            logger.warn(wmsg)
            warnings.warn(wmsg, category=GromacsFailureWarning)

    return {'qtot': qtot,
            'struct': realpath(dirname, 'ionized.gro'),
            'ndx': realpath(dirname, ndx),      # not sure why this is propagated-is it used?
            'mainselection': mainselection,
            }



def solvate(struct='top/protein.pdb', top='top/system.top',
            distance=0.9, boxtype='dodecahedron',
            concentration=0, cation='NA', anion='CL',
            water='tip4p', solvent_name='SOL', with_membrane=False,
            ndx = 'main.ndx', mainselection = '"Protein"',
            dirname='solvate',
            **kwargs):
    """Put protein into box, add water, add counter-ions.

    Currently this really only supports solutes in water. If you need
    to embedd a protein in a membrane then you will require more
    sophisticated approaches.

    However, you *can* supply a protein already inserted in a
    bilayer. In this case you will probably want to set *distance* =
    ``None`` and also enable *with_membrane* = ``True`` (using extra
    big vdw radii for typical lipids).

    .. Note:: The defaults are suitable for solvating a globular
       protein in a fairly tight (increase *distance*!) dodecahedral
       box.

    :Arguments:
      *struct* : filename
          pdb or gro input structure
      *top* : filename
          Gromacs topology
      *distance* : float
          When solvating with water, make the box big enough so that
          at least *distance* nm water are between the solute *struct*
          and the box boundary.
          Set *boxtype*  to ``None`` in order to use a box size in the input
          file (gro or pdb).
      *boxtype* or *bt*: string
          Any of the box types supported by :class:`~gromacs.tools.Editconf`
          (triclinic, cubic, dodecahedron, octahedron). Set the box dimensions
          either with *distance* or the *box* and *angle* keywords.

          If set to ``None`` it will ignore *distance* and use the box
          inside the *struct* file.

          *bt* overrides the value of *boxtype*.
      *box*
          List of three box lengths [A,B,C] that are used by :class:`~gromacs.tools.Editconf`
          in combination with *boxtype* (``bt`` in :program:`editconf`) and *angles*.
          Setting *box* overrides *distance*.
      *angles*
          List of three angles (only necessary for triclinic boxes).
      *concentration* : float
          Concentration of the free ions in mol/l. Note that counter
          ions are added in excess of this concentration.
      *cation* and *anion* : string
          Molecule names of the ions. This depends on the chosen force field.
      *water* : string
          Name of the water model; one of "spc", "spce", "tip3p",
          "tip4p". This should be appropriate for the chosen force
          field. If an alternative solvent is required, simply supply the path to a box
          with solvent molecules (used by :func:`~gromacs.genbox`'s  *cs* argument)
          and also supply the molecule name via *solvent_name*.
      *solvent_name*
          Name of the molecules that make up the solvent (as set in the itp/top).
          Typically needs to be changed when using non-standard/non-water solvents.
          ["SOL"]
      *with_membrane* : bool
           ``True``: use special ``vdwradii.dat`` with 0.1 nm-increased radii on
           lipids. Default is ``False``.
      *ndx* : filename
          How to name the index file that is produced by this function.
      *mainselection* : string
          A string that is fed to :class:`~gromacs.tools.Make_ndx` and
          which should select the solute.
      *dirname* : directory name
          Name of the directory in which all files for the solvation stage are stored.
      *includes*
          List of additional directories to add to the mdp include path
      *kwargs*
          Additional arguments are passed on to
          :class:`~gromacs.tools.Editconf` or are interpreted as parameters to be
          changed in the mdp file.

    """
    sol = solvate_sol(struct=struct, top=top,
                      distance=distance, boxtype=boxtype,
                      water=water, solvent_name=solvent_name, 
                      with_membrane=with_membrane,
                      dirname=dirname, **kwargs)
    
    ion = solvate_ion(struct=sol['struct'], top=top,
                      concentration=concentration, cation=cation, anion=anion,
                      solvent_name=solvent_name, ndx=ndx,
                      mainselection=mainselection, dirname=dirname,
                      **kwargs)
    return ion


def check_mdpargs(d):
    """Check if any arguments remain in dict *d*."""
    if len(d) > 0:
        wmsg = "Unprocessed mdp option are interpreted as options for grompp:\n"+str(d)
        logger.warn(wmsg)
        warnings.warn(wmsg, category=UsageWarning)
    return len(d) == 0

def energy_minimize(dirname='em', mdp=config.templates['em.mdp'],
                    struct='solvate/ionized.gro', top='top/system.top',
                    output='em.pdb', deffnm="em",
                    mdrunner=None, mdrun_args=None,
                    **kwargs):
    """Energy minimize the system.

    This sets up the system (creates run input files) and also runs
    ``mdrun_d``. Thus it can take a while.

    Additional itp files should be in the same directory as the top file.

    Many of the keyword arguments below already have sensible values.

    :Keywords:
       *dirname*
          set up under directory dirname [em]
       *struct*
          input structure (gro, pdb, ...) [solvate/ionized.gro]
       *output*
          output structure (will be put under dirname) [em.pdb]
       *deffnm*
          default name for mdrun-related files [em]
       *top*
          topology file [top/system.top]
       *mdp*
          mdp file (or use the template) [templates/em.mdp]
       *includes*
          additional directories to search for itp files
       *mdrunner*
          :class:`gromacs.run.MDrunner` instance; by default we
          just try :func:`gromacs.mdrun_d` and :func:`gromacs.mdrun` but a
          MDrunner instance gives the user the ability to run mpi jobs
          etc. [None]
       *mdrun_args*
          arguments for *mdrunner* (as a dict), e.g. ``{'nt': 2}``;
          empty by default

          .. versionaddedd:: 0.7.0

       *kwargs*
          remaining key/value pairs that should be changed in the
          template mdp file, eg ``nstxtcout=250, nstfout=250``.

    .. note:: If :func:`~gromacs.mdrun_d` is not found, the function
              falls back to :func:`~gromacs.mdrun` instead.
    """

    structure = realpath(struct)
    topology = realpath(top)
    mdp_template = config.get_template(mdp)
    deffnm = deffnm.strip()

    mdrun_args = {} if mdrun_args is None else mdrun_args

    # write the processed topology to the default output
    kwargs.setdefault('pp', 'processed.top')

    # filter some kwargs that might come through when feeding output
    # from previous stages such as solvate(); necessary because *all*
    # **kwargs must be *either* substitutions in the mdp file *or* valid
    # command line parameters for ``grompp``.
    kwargs.pop('ndx', None)
    # mainselection is not used but only passed through; right now we
    # set it to the default that is being used in all argument lists
    # but that is not pretty. TODO.
    mainselection = kwargs.pop('mainselection', '"Protein"')
    # only interesting when passed from solvate()
    qtot = kwargs.pop('qtot', 0)

    # mdp is now the *output* MDP that will be generated from mdp_template
    mdp = deffnm+'.mdp'
    tpr = deffnm+'.tpr'

    logger.info("[{dirname!s}] Energy minimization of struct={struct!r}, top={top!r}, mdp={mdp!r} ...".format(**vars()))

    cbook.add_mdp_includes(topology, kwargs)

    if qtot != 0:
        # At the moment this is purely user-reported and really only here because
        # it might get fed into the function when using the keyword-expansion pipeline
        # usage paradigm.
        wmsg = "Total charge was reported as qtot = {qtot:g} <> 0; probably a problem.".format(**vars())
        logger.warn(wmsg)
        warnings.warn(wmsg, category=BadParameterWarning)

    with in_dir(dirname):
        unprocessed = cbook.edit_mdp(mdp_template, new_mdp=mdp, **kwargs)
        check_mdpargs(unprocessed)
        gromacs.grompp(f=mdp, o=tpr, c=structure, r=structure, p=topology, **unprocessed)
        mdrun_args.update(v=True, stepout=10, deffnm=deffnm, c=output)
        if mdrunner is None:
            mdrun = run.get_double_or_single_prec_mdrun()
            mdrun(**mdrun_args)
        else:
            if type(mdrunner) is type:
                # class
                # user wants full control and provides simulation.MDrunner **class**
                # NO CHECKING --- in principle user can supply any callback they like
                mdrun = mdrunner(**mdrun_args)
                mdrun.run()
            else:
                # anything with a run() method that takes mdrun arguments...
                try:
                    mdrunner.run(mdrunargs=mdrun_args)
                except AttributeError:
                    logger.error("mdrunner: Provide a gromacs.run.MDrunner class or instance or a callback with a run() method")
                    raise TypeError("mdrunner: Provide a gromacs.run.MDrunner class or instance or a callback with a run() method")

        # em.gro --> gives 'Bad box in file em.gro' warning --- why??
        # --> use em.pdb instead.
        if not os.path.exists(output):
            errmsg = "Energy minimized system NOT produced."
            logger.error(errmsg)
            raise GromacsError(errmsg)
        final_struct = realpath(output)

    logger.info("[{dirname!s}] energy minimized structure {final_struct!r}".format(**vars()))
    return {'struct': final_struct,
            'top': topology,
            'mainselection': mainselection,
            }

def em_schedule(**kwargs):
    """Run multiple energy minimizations one after each other.

    :Keywords:
      *integrators*
           list of integrators (from 'l-bfgs', 'cg', 'steep')
           [['bfgs', 'steep']]
      *nsteps*
           list of maximum number of steps; one for each integrator in
           in the *integrators* list [[100,1000]]
      *kwargs*
           mostly passed to :func:`gromacs.setup.energy_minimize`

    :Returns: dictionary with paths to final structure ('struct') and
              other files

    :Example:
       Conduct three minimizations:
         1. low memory Broyden-Goldfarb-Fletcher-Shannon (BFGS) for 30 steps
         2. steepest descent for 200 steps
         3. finish with BFGS for another 30 steps
       We also do a multi-processor minimization when possible (i.e. for steep
       (and conjugate gradient) by using a :class:`gromacs.run.MDrunner` class
       for a :program:`mdrun` executable compiled for OpenMP in 64 bit (see
       :mod:`gromacs.run` for details)::

          import gromacs.run
          gromacs.setup.em_schedule(struct='solvate/ionized.gro',
                    mdrunner=gromacs.run.MDrunnerOpenMP64,
                    integrators=['l-bfgs', 'steep', 'l-bfgs'],
                    nsteps=[50,200, 50])

    .. Note:: You might have to prepare the mdp file carefully because at the
              moment one can only modify the *nsteps* parameter on a
              per-minimizer basis.
    """

    mdrunner = kwargs.pop('mdrunner', None)
    integrators = kwargs.pop('integrators', ['l-bfgs', 'steep'])
    kwargs.pop('integrator', None)  # clean input; we set intgerator from integrators
    nsteps = kwargs.pop('nsteps', [100, 1000])

    outputs = ['em{0:03d}_{1!s}.pdb'.format(i, integrator) for i,integrator in enumerate(integrators)]
    outputs[-1] = kwargs.pop('output', 'em.pdb')

    files = {'struct': kwargs.pop('struct', None)}  # fake output from energy_minimize()

    for i, integrator in enumerate(integrators):
        struct = files['struct']
        logger.info("[em %d] energy minimize with %s for maximum %d steps", i, integrator, nsteps[i])
        kwargs.update({'struct':struct, 'output':outputs[i],
                       'integrator':integrator, 'nsteps': nsteps[i]})
        if not integrator == 'l-bfgs':
            kwargs['mdrunner'] = mdrunner
        else:
            kwargs['mdrunner'] = None
            logger.warning("[em %d]  Not using mdrunner for L-BFGS because it cannot "
                           "do parallel runs.", i)

        files = energy_minimize(**kwargs)

    return files


def _setup_MD(dirname,
              deffnm='md', mdp=config.templates['md_OPLSAA.mdp'],
              struct=None,
              top='top/system.top', ndx=None,
              mainselection='"Protein"',
              qscript=config.qscript_template, qname=None, startdir=None, mdrun_opts="", budget=None, walltime=1/3.,
              dt=0.002, runtime=1e3, **mdp_kwargs):
    """Generic function to set up a ``mdrun`` MD simulation.

    See the user functions for usage.
    """

    if struct is None:
        raise ValueError('struct must be set to a input structure')
    structure = realpath(struct)
    topology = realpath(top)
    try:
        index = realpath(ndx)
    except AttributeError:  # (that's what realpath(None) throws...)
        index = None        # None is handled fine below

    qname = mdp_kwargs.pop('sgename', qname)    # compatibility for old scripts
    qscript = mdp_kwargs.pop('sge', qscript)    # compatibility for old scripts
    qscript_template = config.get_template(qscript)
    mdp_template = config.get_template(mdp)

    nsteps = int(float(runtime)/float(dt))

    mdp = deffnm + '.mdp'
    tpr = deffnm + '.tpr'
    mainindex = deffnm + '.ndx'
    final_structure = deffnm + '.gro'   # guess... really depends on templates,could also be DEFFNM.pdb

    # write the processed topology to the default output
    mdp_parameters = {'nsteps':nsteps, 'dt':dt, 'pp': 'processed.top'}
    mdp_parameters.update(mdp_kwargs)

    cbook.add_mdp_includes(topology, mdp_parameters)

    logger.info("[%(dirname)s] input mdp  = %(mdp_template)r", vars())
    with in_dir(dirname):
        if not (mdp_parameters.get('Tcoupl','').lower() == 'no' or mainselection is None):
            logger.info("[{dirname!s}] Automatic adjustment of T-coupling groups".format(**vars()))

            # make index file in almost all cases; with mainselection == None the user
            # takes FULL control and also has to provide the template or index
            groups = make_main_index(structure, selection=mainselection,
                                     oldndx=index, ndx=mainindex)
            natoms = {g['name']: float(g['natoms']) for g in groups}
            tc_group_names = ('__main__', '__environment__')   # defined in make_main_index()
            try:
                x = natoms['__main__']/natoms['__environment__']
            except KeyError:
                x = 0   # force using SYSTEM in code below
                wmsg = "Missing __main__ and/or __environment__ index group.\n" \
                       "This probably means that you have an atypical system. You can " \
                       "set mainselection=None and provide your own mdp and index files " \
                       "in order to set up temperature coupling.\n" \
                       "If no T-coupling is required then set Tcoupl='no'.\n" \
                       "For now we will just couple everything to 'System'."
                logger.warn(wmsg)
                warnings.warn(wmsg, category=AutoCorrectionWarning)
            if x < 0.1:
                # couple everything together
                tau_t = firstof(mdp_parameters.pop('tau_t', 0.1))
                ref_t = firstof(mdp_parameters.pop('ref_t', 300))
                # combine all in one T-coupling group
                mdp_parameters['tc-grps'] = 'System'
                mdp_parameters['tau_t'] = tau_t   # this overrides the commandline!
                mdp_parameters['ref_t'] = ref_t   # this overrides the commandline!
                mdp_parameters['gen-temp'] = mdp_parameters.pop('gen_temp', ref_t)
                wmsg = "Size of __main__ is only %.1f%% of __environment__ so " \
                       "we use 'System' for T-coupling and ref_t = %g K and " \
                       "tau_t = %g 1/ps (can be changed in mdp_parameters).\n" \
                       % (x * 100, ref_t, tau_t)
                logger.warn(wmsg)
                warnings.warn(wmsg, category=AutoCorrectionWarning)
            else:
                # couple protein and bath separately
                n_tc_groups = len(tc_group_names)
                tau_t = asiterable(mdp_parameters.pop('tau_t', 0.1))
                ref_t = asiterable(mdp_parameters.pop('ref_t', 300))

                if len(tau_t) != n_tc_groups:
                    tau_t = n_tc_groups * [tau_t[0]]
                    wmsg = "%d coupling constants should have been supplied for tau_t. "\
                        "Using %f 1/ps for all of them." % (n_tc_groups, tau_t[0])
                    logger.warn(wmsg)
                    warnings.warn(wmsg, category=AutoCorrectionWarning)
                if len(ref_t) != n_tc_groups:
                    ref_t = n_tc_groups * [ref_t[0]]
                    wmsg = "%d temperatures should have been supplied for ref_t. "\
                        "Using %g K for all of them." % (n_tc_groups, ref_t[0])
                    logger.warn(wmsg)
                    warnings.warn(wmsg, category=AutoCorrectionWarning)

                mdp_parameters['tc-grps'] = tc_group_names
                mdp_parameters['tau_t'] = tau_t
                mdp_parameters['ref_t'] = ref_t
                mdp_parameters['gen-temp'] = mdp_parameters.pop('gen_temp', ref_t[0])
            index = realpath(mainindex)
        if mdp_parameters.get('Tcoupl','').lower() == 'no':
            logger.info("Tcoupl == no: disabling all temperature coupling mdp options")
            mdp_parameters['tc-grps'] = ""
            mdp_parameters['tau_t'] = ""
            mdp_parameters['ref_t'] = ""
            mdp_parameters['gen-temp'] = ""
        if mdp_parameters.get('Pcoupl','').lower() == 'no':
            logger.info("Pcoupl == no: disabling all pressure coupling mdp options")
            mdp_parameters['tau_p'] = ""
            mdp_parameters['ref_p'] = ""
            mdp_parameters['compressibility'] = ""

        unprocessed = cbook.edit_mdp(mdp_template, new_mdp=mdp, **mdp_parameters)
        check_mdpargs(unprocessed)
        gromacs.grompp(f=mdp, p=topology, c=structure, n=index, o=tpr, **unprocessed)

        runscripts = qsub.generate_submit_scripts(
            qscript_template, deffnm=deffnm, jobname=qname, budget=budget,
            startdir=startdir, mdrun_opts=mdrun_opts, walltime=walltime)

    logger.info("[%(dirname)s] output mdp = %(mdp)r", vars())
    logger.info("[%(dirname)s] output ndx = %(ndx)r", vars())
    logger.info("[%(dirname)s] output tpr = %(tpr)r", vars())
    logger.info("[%(dirname)s] output runscripts = %(runscripts)r", vars())
    logger.info("[%(dirname)s] All files set up for a run time of %(runtime)g ps "
                "(dt=%(dt)g, nsteps=%(nsteps)g)" % vars())

    kwargs = {'struct': realpath(os.path.join(dirname, final_structure)),      # guess
              'top': topology,
              'ndx': index,            # possibly mainindex
              'qscript': runscripts,
              'mainselection': mainselection,
              'deffnm': deffnm,        # return deffnm (tpr = deffnm.tpr!)
              }
    kwargs.update(mdp_kwargs)  # return extra mdp args so that one can use them for prod run
    return kwargs


def MD_restrained(dirname='MD_POSRES', **kwargs):
    """Set up MD with position restraints.

    Additional itp files should be in the same directory as the top file.

    Many of the keyword arguments below already have sensible values. Note that
    setting *mainselection* = ``None`` will disable many of the automated
    choices and is often recommended when using your own mdp file.

    :Keywords:
       *dirname*
          set up under directory dirname [MD_POSRES]
       *struct*
          input structure (gro, pdb, ...) [em/em.pdb]
       *top*
          topology file [top/system.top]
       *mdp*
          mdp file (or use the template) [templates/md.mdp]
       *ndx*
          index file (supply when using a custom mdp)
       *includes*
          additional directories to search for itp files
       *mainselection*
          :program:`make_ndx` selection to select main group ["Protein"]
          (If ``None`` then no canonical index file is generated and
          it is the user's responsibility to set *tc_grps*,
          *tau_t*, and *ref_t* as keyword arguments, or provide the mdp template
          with all parameter pre-set in *mdp* and probably also your own *ndx*
          index file.)
       *deffnm*
          default filename for Gromacs run [md]
       *runtime*
          total length of the simulation in ps [1000]
       *dt*
          integration time step in ps [0.002]
       *qscript*
          script to submit to the queuing system; by default
          uses the template :data:`gromacs.config.qscript_template`, which can
          be manually set to another template from :data:`gromacs.config.templates`;
          can also be a list of template names.
       *qname*
          name to be used for the job in the queuing system [PR_GMX]
       *mdrun_opts*
          option flags for the :program:`mdrun` command in the queuing system
          scripts such as "-stepout 100". [""]
       *kwargs*
          remaining key/value pairs that should be changed in the template mdp
          file, eg ``nstxtcout=250, nstfout=250`` or command line options for
          ``grompp` such as ``maxwarn=1``.

          In particular one can also set **define** and activate
          whichever position restraints have been coded into the itp
          and top file. For instance one could have

             *define* = "-DPOSRES_MainChain -DPOSRES_LIGAND"

          if these preprocessor constructs exist. Note that there
          **must not be any space between "-D" and the value.**

          By default *define* is set to "-DPOSRES".

    :Returns: a dict that can be fed into :func:`gromacs.setup.MD`
              (but check, just in case, especially if you want to
              change the ``define`` parameter in the mdp file)

    .. Note:: The output frequency is drastically reduced for position
              restraint runs by default. Set the corresponding ``nst*``
              variables if you require more output. The `pressure coupling`_
              option *refcoord_scaling* is set to "com" by default (but can
              be changed via *kwargs*) and the pressure coupling
              algorithm itself is set to *Pcoupl* = "Berendsen" to
              run a stable simulation.

    .. _`pressure coupling`: http://manual.gromacs.org/online/mdp_opt.html#pc
    """

    logger.info("[{dirname!s}] Setting up MD with position restraints...".format(**vars()))
    kwargs.setdefault('struct', 'em/em.pdb')
    kwargs.setdefault('qname', 'PR_GMX')
    kwargs.setdefault('define', '-DPOSRES')
    # reduce size of output files
    kwargs.setdefault('nstxout', '50000')   # trr pos
    kwargs.setdefault('nstvout', '50000')   # trr veloc
    kwargs.setdefault('nstfout', '0')       # trr forces
    kwargs.setdefault('nstlog', '500')      # log file
    kwargs.setdefault('nstenergy', '2500')  # edr energy
    kwargs.setdefault('nstxtcout', '5000')  # xtc pos
    # try to get good pressure equilibration
    kwargs.setdefault('refcoord_scaling', 'com')
    kwargs.setdefault('Pcoupl', "Berendsen")

    new_kwargs =  _setup_MD(dirname, **kwargs)

    # clean up output kwargs
    new_kwargs.pop('define', None)          # but make sure that -DPOSRES does not stay...
    new_kwargs.pop('refcoord_scaling', None)
    new_kwargs.pop('Pcoupl', None)
    return new_kwargs

def MD(dirname='MD', **kwargs):
    """Set up equilibrium MD.

    Additional itp files should be in the same directory as the top file.

    Many of the keyword arguments below already have sensible values. Note that
    setting *mainselection* = ``None`` will disable many of the automated
    choices and is often recommended when using your own mdp file.

    :Keywords:
       *dirname*
          set up under directory dirname [MD]
       *struct*
          input structure (gro, pdb, ...) [MD_POSRES/md_posres.pdb]
       *top*
          topology file [top/system.top]
       *mdp*
          mdp file (or use the template) [templates/md.mdp]
       *ndx*
          index file (supply when using a custom mdp)
       *includes*
          additional directories to search for itp files
       *mainselection*
          ``make_ndx`` selection to select main group ["Protein"]
          (If ``None`` then no canonical index file is generated and
          it is the user's responsibility to set *tc_grps*,
          *tau_t*, and *ref_t* as keyword arguments, or provide the mdp template
          with all parameter pre-set in *mdp* and probably also your own *ndx*
          index file.)
       *deffnm*
          default filename for Gromacs run [md]
       *runtime*
          total length of the simulation in ps [1000]
       *dt*
          integration time step in ps [0.002]
       *qscript*
          script to submit to the queuing system; by default
          uses the template :data:`gromacs.config.qscript_template`, which can
          be manually set to another template from :data:`gromacs.config.templates`;
          can also be a list of template names.
       *qname*
          name to be used for the job in the queuing system [MD_GMX]
       *mdrun_opts*
          option flags for the :program:`mdrun` command in the queuing system
          scripts such as "-stepout 100 -dgdl". [""]
       *kwargs*
          remaining key/value pairs that should be changed in the template mdp
          file, e.g. ``nstxtcout=250, nstfout=250`` or command line options for
          :program`grompp` such as ``maxwarn=1``.

    :Returns: a dict that can be fed into :func:`gromacs.setup.MD`
              (but check, just in case, especially if you want to
              change the *define* parameter in the mdp file)
    """

    logger.info("[{dirname!s}] Setting up MD...".format(**vars()))
    kwargs.setdefault('struct', 'MD_POSRES/md.gro')
    kwargs.setdefault('qname', 'MD_GMX')
    return _setup_MD(dirname, **kwargs)



# TODO: autorun (qv MI lipids/setup.pl)

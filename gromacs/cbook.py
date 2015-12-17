# GromacsWrapper -- cbook.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)

"""
:mod:`gromacs.cbook` -- Gromacs Cook Book
=========================================

The :mod:`~gromacs.cbook` (cook book) module contains short recipes for tasks
that are often repeated. In the simplest case this is just one of the
gromacs tools with a certain set of default command line options.

By abstracting and collecting these invocations here, errors can be
reduced and the code snippets can also serve as canonical examples for
how to do simple things.

Miscellaneous canned Gromacs commands
-------------------------------------

Simple commands with new default options so that they solve a specific
problem (see also `Manipulating trajectories and structures`_):

.. function:: rmsd_backbone([s="md.tpr", f="md.xtc"[, ...]])

   Computes the RMSD of the "Backbone" atoms after fitting to the
   "Backbone" (including both translation and rotation).


Manipulating trajectories and structures
----------------------------------------

Standard invocations for manipulating trajectories.

.. function:: trj_compact([s="md.tpr", f="md.xtc", o="compact.xtc"[, ...]])

   Writes an output trajectory or frame with a compact representation
   of the system centered on the protein. It centers on the group
   "Protein" and outputs the whole "System" group.


.. function:: trj_xyfitted([s="md.tpr", f="md.xtc"[, ...]])

   Writes a trajectory centered and fitted to the protein in the XY-plane only.

   This is useful for membrane proteins. The system *must* be oriented so that
   the membrane is in the XY plane. The protein backbone is used for the least
   square fit, centering is done for the whole protein., but this can be
   changed with the *input* = ``('backbone', 'protein','system')`` keyword.

   .. Note:: Gromacs 4.x only

.. autofunction:: trj_fitandcenter
.. autofunction:: cat
.. autoclass:: Frames
   :members:
.. autoclass:: Transformer
   :members:
.. autofunction:: get_volume

Processing output
-----------------

There are cases when a script has to to do different things depending
on the output from a Gromacs tool.

For instance, a common case is to check the total charge after
grompping a tpr file. The ``grompp_qtot`` function does just that.

.. autofunction:: grompp_qtot
.. autofunction:: get_volume
.. autofunction:: parse_ndxlist


Working with topologies and mdp files
-------------------------------------

.. autofunction:: create_portable_topology
.. autofunction:: edit_mdp
.. autofunction:: add_mdp_includes
.. autofunction:: grompp_qtot


Working with index files
------------------------

Manipulation of index files (``ndx``) can be cumbersome because the
``make_ndx`` program is not very sophisticated (yet) compared to
full-fledged atom selection expression as available in Charmm_, VMD_, or
MDAnalysis_. Some tools help in building and interpreting index files.

.. SeeAlso:: The :class:`gromacs.formats.NDX` class can solve a number
             of index problems in a cleaner way than the classes and
             functions here.

.. autoclass:: IndexBuilder
   :members: combine, gmx_resid

.. autofunction:: parse_ndxlist
.. autofunction:: get_ndx_groups
.. autofunction:: make_ndx_captured


.. _MDAnalysis: http://mdanalysis.googlecode.com
.. _VMD: http://www.ks.uiuc.edu/Research/vmd/current/ug/node87.html
.. _Charmm: http://www.charmm.org/html/documentation/c35b1/select.html


File editing functions
----------------------

It is often rather useful to be able to change parts of a template
file. For specialized cases the two following functions are useful:

.. autofunction:: edit_mdp
.. autofunction:: edit_txt
"""

# Right now the simplest thing to do is to just create instances with pre-set
# values; this works fine and is succinct but has some disadvantages:
# * the underlying gromacs tool is executed to extract the help string; this
#   adds to the import time
# * adding documentation is awkward
#
# For more complicated cases one is probably better off by properly deriving a
# new class and set arguments explicitly in init (using kwargs['flag'] =
# default) ... or I can write some meta(??) class to do this nicely

from __future__ import absolute_import, with_statement

__docformat__ = "restructuredtext en"

import sys
import os
import re
import warnings
import tempfile
import shutil
import glob

import logging
logger = logging.getLogger('gromacs.cbook')

import gromacs
from .exceptions import GromacsError, BadParameterWarning, MissingDataWarning, GromacsValueWarning, GromacsImportWarning
from . import tools
from . import utilities
from .utilities import asiterable

def _define_canned_commands():
    """Define functions for the top level name space.

    Definitions are collected here so that they can all be wrapped in
    a try-except block that avoids code failing when the Gromacs tools
    are not available --- in some cases they are not necessary to use
    parts of GromacsWrapper.

    .. Note:: Any function defined here **must be listed in ``global``**!
    """
    global trj_compact, rmsd_backbone, trj_fitted, trj_xyfitted

    trj_compact = tools.Trjconv(ur='compact', center=True, boxcenter='tric', pbc='mol',
                                input=('protein','system'),
                                doc="""
Writes a compact representation of the system centered on the protein""")

    rmsd_backbone = tools.G_rms(what='rmsd', fit='rot+trans',
                                input=('Backbone','Backbone'),
                                doc="""
Computes RMSD of backbone after fitting to the backbone.""")

    trj_fitted = tools.Trjconv(fit='rot+trans',
                               input=('backbone', 'system'),
                               doc="""
Writes a trajectory fitted to the protein backbone.

Note that this does *not* center; if center is required, the *input*
selection should have the group to be centered on in second position,
e.g. ``input = ('backbone', 'Protein', System')``.
""")


    # Gromacs 4.x
    trj_xyfitted = tools.Trjconv(fit='rotxy+transxy',
                                 input=('backbone', 'protein','system'),
                                 doc="""
Writes a trajectory fitted to the protein in the XY-plane only.

This is useful for membrane proteins. The system *must* be oriented so
that the membrane is in the XY plane. The protein backbone is used
for the least square fit, centering is done for the whole protein.

Note that centering together with fitting does not always work well
and that one sometimes need two runs of trjconv: one to center and
one to fit.

.. Note:: Gromacs 4.x only""")

    # end of _define_canned_commands

try:
    _define_canned_commands()
except (OSError, ImportError, GromacsError):
    warnings.warn("Failed to define a number of commands in gromacs.cbook. Most likely the "
                  "Gromacs installation cannot be found --- source GMXRC!",
                  category=GromacsImportWarning)
    logger.error("Failed to define a number of commands in gromacs.cbook. Most likely the "
                  "Gromacs installation cannot be found --- source GMXRC!")
finally:
    del _define_canned_commands

def trj_fitandcenter(xy=False, **kwargs):
    """Center everything and make a compact representation (pass 1) and fit the system to a reference (pass 2).

    :Keywords:
       *s*
           input structure file (tpr file required to make molecule whole);
           if a list or tuple is provided then s[0] is used for pass 1 (should be a tpr)
           and s[1] is used for the fitting step (can be a pdb of the whole system)

           If a second structure is supplied then it is assumed that the fitted
           trajectory should *not* be centered.
       *f*
           input trajectory
       *o*
           output trajectory
       *input*
           A list with three groups. The default is
               ['backbone', 'protein','system']
           The fit command uses all three (1st for least square fit,
           2nd for centering, 3rd for output), the centered/make-whole stage use
           2nd for centering and 3rd for output.
       *input1*
           If *input1* is supplied then *input* is used exclusively
           for the fitting stage (pass 2) and *input1* for the centering (pass 1).
       *n*
           Index file used for pass 1 and pass 2.
       *n1*
           If *n1* is supplied then index *n1* is only used for pass 1
           (centering) and *n* for pass 2 (fitting).
       *xy* : boolean
           If ``True`` then only do a rot+trans fit in the xy plane
           (good for membrane simulations); default is ``False``.
       *kwargs*
           All other arguments are passed to :class:`~gromacs.tools.Trjconv`.

    Note that here we first center the protein and create a compact box, using
    ``-pbc mol -ur compact -center -boxcenter tric`` and write an intermediate
    xtc. Then in a second pass we perform a rotation+translation fit (or
    restricted to the xy plane if *xy* = ``True`` is set) on the intermediate
    xtc to produce the final trajectory. Doing it in this order has the
    disadvantage that the solvent box is rotating around the protein but the
    opposite order (with center/compact second) produces strange artifacts
    where columns of solvent appear cut out from the box---it probably means
    that after rotation the information for the periodic boundaries is not
    correct any more.

    Most kwargs are passed to both invocations of
    :class:`gromacs.tools.Trjconv` so it does not really make sense to use eg
    *skip*; in this case do things manually.

    By default the *input* to the fit command is ('backbone',
    'protein','system'); the compact command always uses the second and third
    group for its purposes or if this fails, prompts the user.

    Both steps cannot performed in one pass; this is a known limitation of
    ``trjconv``. An intermediate temporary XTC files is generated which should
    be automatically cleaned up unless bad things happened.

    The function tries to honour the input/output formats. For instance, if you
    want trr output you need to supply a trr file as input and explicitly give
    the output file also a trr suffix.

    .. Note:: For big trajectories it can **take a very long time**
              and consume a **large amount of temporary diskspace**.

    We follow the `g_spatial documentation`_ in preparing the trajectories::

       trjconv -s a.tpr -f a.xtc -o b.xtc -center -boxcenter tric -ur compact -pbc mol
       trjconv -s a.tpr -f b.xtc -o c.xtc -fit rot+trans

    .. _`g_spatial documentation`: http://www.gromacs.org/Documentation/Gromacs_Utilities/g_spatial
    """
    if xy:
        fitmode = 'rotxy+transxy'
        kwargs.pop('fit', None)
    else:
        fitmode = kwargs.pop('fit', 'rot+trans')  # user can use progressive, too

    intrj = kwargs.pop('f', None)
    # get the correct suffix for the intermediate step: only trr will
    # keep velocities/forces!
    suffix = os.path.splitext(intrj)[1]
    if not suffix in ('xtc', 'trr'):
        suffix = '.xtc'
    outtrj = kwargs.pop('o', None)

    ndx = kwargs.pop('n', None)
    ndxcompact = kwargs.pop('n1', ndx)

    structures = kwargs.pop('s', None)
    if type(structures) in (tuple, list):
        try:
            compact_structure, fit_structure = structures
        except:
            raise ValueError("argument s must be a pair of tpr/pdb files or a single structure file")
    else:
        compact_structure = fit_structure = structures


    inpfit = kwargs.pop('input', ('backbone', 'protein','system'))
    try:
        _inpcompact = inpfit[1:]     # use 2nd and 3rd group for compact
    except TypeError:
        _inpcompact = None
    inpcompact = kwargs.pop('input1', _inpcompact)  # ... or the user supplied ones

    fd, tmptrj = tempfile.mkstemp(suffix=suffix, prefix='pbc_compact_')

    logger.info("Input structure for PBC:  %(compact_structure)r" % vars())
    logger.info("Input structure for fit:  %(fit_structure)r" % vars())
    logger.info("Input trajectory:  %(intrj)r" % vars())
    logger.info("Output trajectory: %(outtrj)r"% vars())
    logger.debug("Writing temporary trajectory %(tmptrj)r (will be auto-cleaned)." % vars())
    sys.stdout.flush()
    try:
        gromacs.trjconv(s=compact_structure, f=intrj, o=tmptrj, n=ndxcompact,
                        ur='compact', center=True, boxcenter='tric', pbc='mol',
                        input=inpcompact, **kwargs)
        # explicitly set pbc="none" for the fitting stage (anything else will produce rubbish and/or
        # complaints from Gromacs)
        kwargs['pbc'] = "none"
        if compact_structure == fit_structure:
            # fit as ususal, including centering
            # (Is center=True really necessary? -- note, if I remove center=True then
            # I MUST fiddle inpfit as below!!)
            gromacs.trjconv(s=fit_structure, f=tmptrj, o=outtrj, n=ndx, fit=fitmode, center=True, input=inpfit, **kwargs)
        else:
            # make sure that we fit EXACTLY as the user wants
            inpfit = [inpfit[0], inpfit[-1]]
            gromacs.trjconv(s=fit_structure, f=tmptrj, o=outtrj, n=ndx, fit=fitmode, input=inpfit, **kwargs)
    finally:
        utilities.unlink_gmx(tmptrj)

def cat(prefix="md", dirname=os.path.curdir, partsdir="parts", fulldir="full",
        resolve_multi="pass"):
    """Concatenate all parts of a simulation.

    The xtc, trr, and edr files in *dirname* such as prefix.xtc,
    prefix.part0002.xtc, prefix.part0003.xtc, ... are

       1) moved to the *partsdir* (under *dirname*)
       2) concatenated with the Gromacs tools to yield prefix.xtc, prefix.trr,
          prefix.edr, prefix.gro (or prefix.md) in *dirname*
       3) Store these trajectories in *fulldir*

    .. Note:: Trajectory files are *never* deleted by this function to avoid
              data loss in case of bugs. You will have to clean up yourself
              by deleting *dirname*/*partsdir*.

              Symlinks for the trajectories are *not* handled well and
              break the function. Use hard links instead.

    .. Warning:: If an exception occurs when running this function then make
                 doubly and triply sure where your files are before running
                 this function again; otherwise you might **overwrite data**.
                 Possibly you will need to manually move the files from *partsdir*
                 back into the working directory *dirname*; this should onlu overwrite
                 generated files so far but *check carefully*!

    :Keywords:
        *prefix*
            deffnm of the trajectories [md]
        *resolve_multi"
            how to deal with multiple "final" gro or pdb files: normally there should
            only be one but in case of restarting from the checkpoint of a finished
            simulation one can end up with multiple identical ones.
              - "pass" : do nothing and log a warning
              - "guess" : take prefix.pdb or prefix.gro if it exists, otherwise the one of
                          prefix.partNNNN.gro|pdb with the highes NNNN
        *dirname*
            change to *dirname* and assume all tarjectories are located there [.]
        *partsdir*
             directory where to store the input files (they are moved out of the way);
             *partsdir* must be manually deleted [parts]
        *fulldir*
             directory where to store the final results [full]
    """

    gmxcat = {'xtc': gromacs.trjcat,
              'trr': gromacs.trjcat,
              'edr': gromacs.eneconv,
              'log': utilities.cat,
              }

    def _cat(prefix, ext, partsdir=partsdir, fulldir=fulldir):
        filenames = glob_parts(prefix, ext)
        if ext.startswith('.'):
            ext = ext[1:]
        outfile = os.path.join(fulldir, prefix + '.' + ext)
        if not filenames:
            return None
        nonempty_files = []
        for f in filenames:
            if os.stat(f).st_size == 0:
                logger.warn("File %(f)r is empty, skipping" % vars())
                continue
            if os.path.islink(f):
                # TODO: re-write the symlink to point to the original file
                errmsg = "Symbolic links do not work (file %(f)r), sorry. " \
                    "CHECK LOCATION OF FILES MANUALLY BEFORE RUNNING gromacs.cbook.cat() AGAIN!" % vars()
                logger.exception(errmsg)
                raise NotImplementedError(errmsg)
            shutil.move(f, partsdir)
            nonempty_files.append(f)
        filepaths = [os.path.join(partsdir, f) for f in nonempty_files]
        gmxcat[ext](f=filepaths, o=outfile)
        return outfile

    _resolve_options = ("pass", "guess")
    if not resolve_multi in _resolve_options:
        raise ValueError("resolve_multi must be one of %(_resolve_options)r, "
                         "not %(resolve_multi)r" % vars())

    if fulldir == os.path.curdir:
        wmsg = "Using the current directory as fulldir can potentially lead to data loss if you run this function multiple times."
        logger.warning(wmsg)
        warnings.warn(wmsg, category=BadParameterWarning)

    with utilities.in_dir(dirname, create=False):
        utilities.mkdir_p(partsdir)
        utilities.mkdir_p(fulldir)
        for ext in ('log', 'edr', 'trr', 'xtc'):
            logger.info("[%(dirname)s] concatenating %(ext)s files...", vars())
            outfile = _cat(prefix, ext, partsdir)
            logger.info("[%(dirname)s] created %(outfile)r", vars())
        for ext in ('gro', 'pdb'):              # XXX: ugly, make method out of parts?
            filenames = glob_parts(prefix, ext)
            if len(filenames) == 0:
                continue                        # goto next ext
            elif len(filenames) == 1:
                pick = filenames[0]
            else:
                if resolve_multi == "pass":
                    logger.warning("[%(dirname)s] too many output structures %(filenames)r, "
                                   "cannot decide which one --- resolve manually!", vars())
                    for f in filenames:
                        shutil.move(f, partsdir)
                    continue                    # goto next ext
                elif resolve_multi == "guess":
                    pick = prefix + '.' + ext
                    if not pick in filenames:
                        pick = filenames[-1]  # filenames are ordered with highest parts at end
            final = os.path.join(fulldir, prefix + '.' + ext)
            shutil.copy(pick, final)  # copy2 fails on nfs with Darwin at least
            for f in filenames:
                shutil.move(f, partsdir)
            logger.info("[%(dirname)s] collected final structure %(final)r "
                        "(from %(pick)r)", vars())


    partsdirpath = utilities.realpath(dirname, partsdir)
    logger.warn("[%(dirname)s] cat() complete in %(fulldir)r but original files "
                "in %(partsdirpath)r must be manually removed", vars())

def glob_parts(prefix, ext):
    """Find files from a continuation run"""
    if ext.startswith('.'):
        ext = ext[1:]
    files = glob.glob(prefix+'.'+ext) + glob.glob(prefix+'.part[0-9][0-9][0-9][0-9].'+ext)
    files.sort()   # at least some rough sorting...
    return files


class Frames(object):
    """A iterator that transparently provides frames from a trajectory.

    The iterator chops a trajectory into individual frames for
    analysis tools that only work on separate structures such as
    ``gro`` or ``pdb`` files. Instead of turning the whole trajectory
    immediately into pdb files (and potentially filling the disk), the
    iterator can be instructed to only provide a fixed number of
    frames and compute more frames when needed.

    .. Note:: Setting a limit on the number of frames on disk can lead
              to longish waiting times because ``trjconv`` must
              re-seek to the middle of the trajectory and the only way
              it can do this at the moment is by reading frames
              sequentially. This might still be preferrable to filling
              up a disk, though.

    .. Warning:: The *maxframes* option is not implemented yet; use
                 the *dt* option or similar to keep the number of frames
                 manageable.
    """

    def __init__(self, structure, trj, maxframes=None, format='pdb', **kwargs):
        """Set up the Frames iterator.

        :Arguments:
           structure
              name of a structure file (tpr, pdb, ...)
           trj
              name of the trajectory (xtc, trr, ...)
           format
              output format for the frames, eg "pdb" or "gro" [pdb]
           maxframes : int
             maximum number of frames that are extracted to disk at
             one time; set to ``None`` to extract the whole trajectory
             at once. [``None``]
           kwargs
             All other arguments are passed to
             `class:~gromacs.tools.Trjconv`; the only options that
             cannot be changed are *sep* and the output file name *o*.

        """
        self.structure = structure  # tpr or equivalent
        self.trj = trj              # xtc, trr, ...
        self.maxframes = maxframes
        if not self.maxframes is None:
            raise NotImplementedError('sorry, maxframes feature not implemented yet')

        self.framedir = tempfile.mkdtemp(prefix="Frames_", suffix='_'+format)
        self.frameprefix = os.path.join(self.framedir, 'frame')
        self.frametemplate = self.frameprefix + '%d' + '.' + format  # depends on trjconv
        self.frameglob = self.frameprefix + '*' + '.' + format
        kwargs['sep'] = True
        kwargs['o'] = self.frameprefix + '.' + format
        kwargs.setdefault('input', ('System',))
        self.extractor = tools.Trjconv(s=self.structure, f=self.trj, **kwargs)

        #: Holds the current frame number of the currently extracted
        #: batch of frames. Increases when iterating.
        self.framenumber = 0
        #: Total number of frames read so far; only important when *maxframes* > 0 is used.
        self.totalframes = 0

    def extract(self):
        """Extract frames from the trajectory to the temporary directory."""
        # XXX: extract everything at the moment, logic for maxframes not done yet
        self.extractor.run()

    @property
    def all_frames(self):
        """Unordered list of all frames currently held on disk."""
        return glob.glob(self.frameglob)

    @property
    def current_framename(self):
        return self.frametemplate % self.framenumber

    def __iter__(self):
        """Primitive iterator."""
        frames = self.all_frames
        if len(frames) == 0:
            self.extract()
            frames = self.all_frames

        # filenames are 'Frame0.pdb', 'Frame11.pdb', ... so I must
        # order manually because glob does not give it in sequence.
        for i in xrange(len(frames)):
            self.framenumber = i
            yield self.current_framename
        self.totalframes += len(frames)

    def delete_frames(self):
        """Delete all frames."""
        for frame in glob.glob(self.frameglob):
            os.unlink(frame)

    def cleanup(self):
        """Clean up all temporary frames (which can be HUGE)."""
        shutil.rmtree(self.framedir)
        self.framedir = None

    def __del__(self):
        if not self.framedir is None:
            self.cleanup()

# Working with topologies
# -----------------------

# grompp that does not raise an exception; setting up runs the command to get the docs so
# we only want to do this once at the module level and not inside a function that can be called
# repeatedly
grompp_warnonly = tools.Grompp(failure="warn",
                               doc="grompp wrapper that only warns on failure but does not raise :exc:`GromacsError`")

def grompp_qtot(*args, **kwargs):
    """Run ``gromacs.grompp`` and return the total charge of the  system.

    :Arguments:
       The arguments are the ones one would pass to :func:`gromacs.grompp`.
    :Returns:
       The total charge as reported

    Some things to keep in mind:

    * The stdout output of grompp is only shown when an error occurs. For
      debugging, look at the log file or screen output and try running the
      normal :func:`gromacs.grompp` command and analyze the output if the
      debugging messages are not sufficient.

    * Check that ``qtot`` is correct. Because the function is based on pattern
      matching of the informative output of :program:`grompp` it can break when
      the output format changes. This version recognizes lines like ::

            '  System has non-zero total charge: -4.000001e+00'

      using the regular expression
      :regexp:`System has non-zero total charge: *(?P<qtot>[-+]?\d*\.\d+([eE][-+]\d+)?)`.

    """
    qtot_pattern = re.compile('System has non-zero total charge: *(?P<qtot>[-+]?\d*\.\d+([eE][-+]\d+)?)')
    # make sure to capture ALL output
    kwargs['stdout'] = False
    kwargs['stderr'] = False
    rc, output, error = grompp_warnonly(*args, **kwargs)
    gmxoutput = "\n".join([x for x in [output, error] if x is not None])
    if rc != 0:
        # error occured and we want to see the whole output for debugging
        msg = "grompp_qtot() failed. See warning and screen output for clues."
        logger.error(msg)
        import sys
        sys.stderr.write("=========== grompp (stdout/stderr) ============\n")
        sys.stderr.write(gmxoutput)
        sys.stderr.write("===============================================\n")
        sys.stderr.flush()
        raise GromacsError(rc, msg)
    qtot = 0
    for line in gmxoutput.split('\n'):
        m = qtot_pattern.search(line)
        if m:
            qtot = float(m.group('qtot'))
            break
    logger.info("system total charge qtot = %(qtot)r" % vars())
    return qtot

def _mdp_include_string(dirs):
    """Generate a string that can be added to a mdp 'include = ' line."""
    include_paths = [os.path.expanduser(p) for p in dirs]
    return ' -I'.join([''] + include_paths)

def add_mdp_includes(topology=None, kwargs=None):
    """Set the mdp *include* key in the *kwargs* dict.

    1. Add the directory containing *topology*.
    2. Add all directories appearing under the key *includes*
    3. Generate a string of the form "-Idir1 -Idir2 ..." that
       is stored under the key *include* (the corresponding
       mdp parameter)

    By default, the directories ``.`` and ``..`` are also added to the
    *include* string for the mdp; when fed into
    :func:`gromacs.cbook.edit_mdp` it will result in a line such as ::

      include = -I. -I.. -I../topology_dir ....

    Note that the user can always override the behaviour by setting
    the *include* keyword herself; in this case this function does
    nothing.

    If no *kwargs* were supplied then a dict is generated with the
    single *include* entry.

    :Arguments:
       *topology* : top filename
          Topology file; the name of the enclosing directory is added
          to the include path (if supplied) [``None``]
       *kwargs* : dict
          Optional dictionary of mdp keywords; will be modified in place.
          If it contains the *includes* keyword with either a single string
          or a list of strings then these paths will be added to the
          include statement.
    :Returns: *kwargs* with the *include* keyword added if it did not
              exist previously; if the keyword already existed, nothing
              happens.

    .. Note:: The *kwargs* dict is **modified in place**. This
              function is a bit of a hack. It might be removed once
              all setup functions become methods in a nice class.
    """
    if kwargs is None:
        kwargs = {}

    include_dirs = ['.', '..']      # should . & .. always be added?
    if topology is not None:
        # half-hack: find additional itps in the same directory as the topology
        topology_dir = os.path.dirname(topology)
        include_dirs.append(topology_dir)

    include_dirs.extend(asiterable(kwargs.pop('includes', [])))  # includes can be a list or a string

    # 1. setdefault: we do nothing if user defined include
    # 2. modify input in place!
    kwargs.setdefault('include', _mdp_include_string(include_dirs))
    return kwargs

def filter_grompp_options(**kwargs):
    """Returns one dictionary only containing valid :program:`grompp` options and everything else.

    Option list is hard coded and nased on :class:`~gromacs.tools.grompp` 4.5.3.

    :Returns: ``(grompp_dict, other_dict)``

    .. versionadded:: 0.2.4
    """
    grompp_options = ('f','po','c','r','rb','n','p','pp','o','t','e',  # files
                      'h', 'noh', 'version', 'noversion', 'nice', 'v', 'nov',
                      'time', 'rmvsbds', 'normvsbds', 'maxwarn', 'zero', 'nozero',
                      'renum', 'norenum')
    grompp = dict((k,v) for k,v in kwargs.items() if k in grompp_options)
    other =  dict((k,v) for k,v in kwargs.items() if k not in grompp_options)
    return grompp, other

def create_portable_topology(topol, struct, **kwargs):
    """Create a processed topology.

    The processed (or portable) topology file does not contain any
    ``#include`` statements and hence can be easily copied around. It
    also makes it possible to re-grompp without having any special itp
    files available.

    :Arguments:
      *topol*
          topology file
      *struct*
          coordinat (structure) file

    :Keywords:
      *processed*
          name of the new topology file; if not set then it is named like
          *topol* but with ``pp_`` prepended
      *includes*
          path or list of paths of directories in which itp files are
          searched for
      *grompp_kwargs**
          other options for :program:`grompp` such as ``maxwarn=2`` can
          also be supplied

    :Returns: full path to the processed topology
    """
    _topoldir, _topol = os.path.split(topol)
    processed = kwargs.pop('processed', os.path.join(_topoldir, 'pp_'+_topol))
    grompp_kwargs, mdp_kwargs = filter_grompp_options(**kwargs)
    mdp_kwargs = add_mdp_includes(topol, mdp_kwargs)
    with tempfile.NamedTemporaryFile(suffix='.mdp') as mdp:
        mdp.write('; empty mdp file\ninclude = %(include)s\n' % mdp_kwargs)
        mdp.flush()
        grompp_kwargs['p'] = topol
        grompp_kwargs['pp'] = processed
        grompp_kwargs['f'] =  mdp.name
        grompp_kwargs['c'] = struct
        grompp_kwargs['v'] = False
        try:
            gromacs.grompp(**grompp_kwargs)
        finally:
            utilities.unlink_gmx('topol.tpr', 'mdout.mdp')
    return utilities.realpath(processed)

def get_volume(f):
    """Return the volume in nm^3 of structure file *f*.

    (Uses :func:`gromacs.editconf`; error handling is not good)
    """
    fd, temp = tempfile.mkstemp('.gro')
    try:
        rc,out,err = gromacs.editconf(f=f, o=temp, stdout=False)
    finally:
        os.unlink(temp)
    return [float(x.split()[1]) for x in out.splitlines()
            if x.startswith('Volume:')][0]


# Editing textual input files
# ---------------------------

def edit_mdp(mdp, new_mdp=None, extend_parameters=None, **substitutions):
    """Change values in a Gromacs mdp file.

    Parameters and values are supplied as substitutions, eg ``nsteps=1000``.

    By default the template mdp file is **overwritten in place**.

    If a parameter does not exist in the template then it cannot be substituted
    and the parameter/value pair is returned. The user has to check the
    returned list in order to make sure that everything worked as expected. At
    the moment it is not possible to automatically append the new values to the
    mdp file because of ambiguities when having to replace dashes in parameter
    names with underscores (see the notes below on dashes/underscores).

    If a parameter is set to the value ``None`` then it will be ignored.

    :Arguments:
        *mdp* : filename
            filename of input (and output filename of ``new_mdp=None``)
        *new_mdp* : filename
            filename of alternative output mdp file [None]
        *extend_parameters* : string or list of strings
            single parameter or list of parameters for which the new values
            should be appended to the existing value in the mdp file. This
            makes mostly sense for a single parameter, namely 'include', which
            is set as the default. Set to ``[]`` to disable. ['include']
        *substitutions*
            parameter=value pairs, where parameter is defined by the Gromacs
            mdp file; dashes in parameter names have to be replaced by
            underscores. If a value is a list-like object then the items are
            written as a sequence, joined with spaces, e.g. ::

               ref_t=[310,310,310] --->  ref_t = 310 310 310

    :Returns:
        Dict of parameters that have *not* been substituted.

    **Example** ::

       edit_mdp('md.mdp', new_mdp='long_md.mdp', nsteps=100000, nstxtcout=1000, lincs_iter=2)

    .. Note::

       * Dashes in Gromacs mdp parameters have to be replaced by an underscore
         when supplied as python keyword arguments (a limitation of python). For example
         the MDP syntax is  ``lincs-iter = 4`` but the corresponding  keyword would be
         ``lincs_iter = 4``.
       * If the keyword is set as a dict key, eg ``mdp_params['lincs-iter']=4`` then one
         does not have to substitute.
       * Parameters *aa_bb* and *aa-bb* are considered the same (although this should
         not be a problem in practice because there are no mdp parameters that only
         differ by a underscore).
       * This code is more compact in ``Perl`` as one can use ``s///`` operators:
         ``s/^(\s*${key}\s*=\s*).*/$1${val}/``

    .. SeeAlso:: One can also load the mdp file with
                :class:`gromacs.formats.MDP`, edit the object (a dict), and save it again.
    """
    if new_mdp is None:
        new_mdp = mdp
    if extend_parameters is None:
        extend_parameters = ['include']
    else:
        extend_parameters = list(asiterable(extend_parameters))

    # None parameters should be ignored (simple way to keep the template defaults)
    substitutions = {k: v for k,v in substitutions.items() if not v is None}

    params = substitutions.keys()[:]   # list will be reduced for each match

    def demangled(p):
        """Return a RE string that matches the parameter."""
        return p.replace('_', '[-_]')  # must catch either - or _

    patterns = {parameter:
                      re.compile("""\
                       (?P<assignment>\s*%s\s*=\s*)  # parameter == everything before the value
                       (?P<value>[^;]*)              # value (stop before comment=;)
                       (?P<comment>\s*;.*)?          # optional comment
                       """ % demangled(parameter), re.VERBOSE)
                     for parameter in substitutions}

    target = tempfile.TemporaryFile()
    with open(mdp) as src:
        logger.info("editing mdp = %r: %r" % (mdp, substitutions.keys()))
        for line in src:
            new_line = line.strip()  # \n must be stripped to ensure that new line is built without break
            for p in params[:]:
                m = patterns[p].match(new_line)
                if m:
                    # I am too stupid to replace a specific region in the string so I rebuild it
                    # (matching a line and then replacing value requires TWO re calls)
                    #print 'line:' + new_line
                    #print m.groupdict()
                    if m.group('comment') is None:
                        comment = ''
                    else:
                        comment = " "+m.group('comment')
                    assignment = m.group('assignment')
                    if not assignment.endswith(' '):
                        assignment += ' '
                    # build new line piece-wise:
                    new_line = assignment
                    if p in extend_parameters:
                        # keep original value and add new stuff at end
                        new_line += str(m.group('value')) + ' '
                    # automatically transform lists into space-separated string values
                    value = " ".join(map(str, asiterable(substitutions[p])))
                    new_line += value + comment
                    params.remove(p)
                    break
            target.write(new_line+'\n')
    target.seek(0)
    # XXX: Is there a danger of corrupting the original mdp if something went wrong?
    with open(new_mdp, 'w') as final:
        shutil.copyfileobj(target, final)
    target.close()
     # return all parameters that have NOT been substituted
    if len(params) > 0:
        logger.warn("Not substituted in %(new_mdp)r: %(params)r" % vars())
    return {p: substitutions[p] for p in params}

def edit_txt(filename, substitutions, newname=None):
    """Primitive text file stream editor.

    This function can be used to edit free-form text files such as the
    topology file. By default it does an **in-place edit** of
    *filename*. If *newname* is supplied then the edited
    file is written to *newname*.

    :Arguments:
       *filename*
           input text file
       *substitutions*
           substitution commands (see below for format)
       *newname*
           output filename; if ``None`` then *filename* is changed in
           place [``None``]

    *substitutions* is a list of triplets; the first two elements are regular
    expression strings, the last is the substitution value. It mimics
    ``sed`` search and replace. The rules for *substitutions*:

    .. productionlist::
       substitutions: "[" search_replace_tuple, ... "]"
       search_replace_tuple: "(" line_match_RE "," search_RE "," replacement ")"
       line_match_RE: regular expression that selects the line (uses match)
       search_RE: regular expression that is searched in the line
       replacement: replacement string for search_RE

    Running :func:`edit_txt` does pretty much what a simple ::

         sed /line_match_RE/s/search_RE/replacement/

    with repeated substitution commands does.

    Special replacement values:
    - ``None``: the rule is ignored
    - ``False``: the line is deleted (even if other rules match)

    .. note::

       * No sanity checks are performed and the substitutions must be supplied
         exactly as shown.
       * All substitutions are applied to a line; thus the order of the substitution
         commands may matter when one substitution generates a match for a subsequent rule.
       * If replacement is set to ``None`` then the whole expression is ignored and
         whatever is in the template is used. To unset values you must provided an
         empty string or similar.
       * Delete a matching line if replacement=``False``.
    """
    if newname is None:
        newname = filename

    # No sanity checks (figure out later how to give decent diagnostics).
    # Filter out any rules that have None in replacement.
    _substitutions = [{'lRE': re.compile(str(lRE)),
                       'sRE': re.compile(str(sRE)),
                       'repl': repl}
                      for lRE,sRE,repl in substitutions if not repl is None]

    target = tempfile.TemporaryFile()
    with open(filename) as src:
        logger.info("editing txt = %r (%d substitutions)" % (filename, len(substitutions)))
        for line in src:
            keep_line = True
            for subst in _substitutions:
                m = subst['lRE'].match(line)
                if m:              # apply substition to this line?
                    logger.debug('match:    '+line.rstrip())
                    if subst['repl'] is False:   # special rule: delete line
                        keep_line = False
                    else:                        # standard replacement
                        line = subst['sRE'].sub(str(subst['repl']), line)
                        logger.debug('replaced: '+line.rstrip())
            if keep_line:
                target.write(line)
            else:
                logger.debug("Deleting line %r", line)

    target.seek(0L)
    with open(newname, 'w') as final:
        shutil.copyfileobj(target, final)
    target.close()
    logger.info("edited txt = %(newname)r" % vars())

def remove_molecules_from_topology(filename, **kwargs):
    """Remove autogenerated [ molecules ] entries from *filename*.

    Valid entries in ``[ molecules ]`` below the default *marker*
    are removed. For example, a topology file such as ::

       [ molecules ]
       Protein    1
       SOL      213
       ; The next line is the marker!
       ; Gromacs auto-generated entries follow:
       SOL            12345
       NA+     15
       CL-      16
       ; This is a comment that is NOT deleted.
       SOL            333

    would become::

       [ molecules ]
       Protein    1
       SOL      213
       ; The next line is the marker!
       ; Gromacs auto-generated entries follow:
       ; This is a comment that is NOT deleted.

    Valid molecule lines look like ``SOL 1234``, ``NA 17`` etc. The
    actual regular expression used is "\s*[\w+_-]+\s+\d+\s*(;.*)?$".

    In order to use this function, the marker line has to be manually
    added to the topology file.

    :Arguments:
      *filename*
         The topology file that includes the  ``[ molecules ]`` section.
         It is **edited in place**.
      *marker*
         Any ``[ molecules ]`` entries below this pattern (python regular
         expression) are removed. Leading white space is ignored. ``None``
         uses the default as described above.
    """
    marker = kwargs.pop('marker', None)
    if marker is None:
        marker = "; Gromacs auto-generated entries follow:"
    logger.debug("Scrubbed [ molecules ]: marker = %(marker)r", vars())

    p_marker = re.compile("\s*%s" % marker)
    p_molecule = re.compile("\s*[\w+_-]+\s+\d+\s*(;.*)?$")
    target = tempfile.TemporaryFile()
    with open(filename) as src:
        autogenerated = False
        n_removed = 0
        for line in src:
            if p_marker.match(line):
                autogenerated = True
            if autogenerated and p_molecule.match(line):
                n_removed += 1
                continue  # remove by skipping
            target.write(line)
    if autogenerated and n_removed > 0:
        target.seek(0L)
        with open(filename, 'w') as final:   # overwrite original!
            shutil.copyfileobj(target, final)
        logger.info("Removed %(n_removed)d autogenerated [ molecules ] from "
                    "topol = %(filename)r" % vars())
    target.close()
    return n_removed


# Working with index files and index groups
# -----------------------------------------
#

#: compiled regular expression to match a list of index groups
#: in the output of ``make_ndx``s <Enter> (empty) command.
NDXLIST = re.compile(r""">\s+\n    # '> ' marker line from '' input (input not echoed)
                 \n                # empty line
                 (?P<LIST>         # list of groups
                  (                # consists of repeats of the same pattern:
                    \s*\d+         # group number
                    \s+[^\s]+\s*:  # group name, separator ':'
                    \s*\d+\satoms  # number of atoms in group
                    \n
                   )+              # multiple repeats
                  )""", re.VERBOSE)
#: compiled regular expression to match a single line of
#: ``make_ndx`` output (e.g. after a successful group creation)
NDXGROUP = re.compile(r"""
                     \s*(?P<GROUPNUMBER>\d+)      # group number
                     \s+(?P<GROUPNAME>[^\s]+)\s*: # group name, separator ':'
                     \s*(?P<NATOMS>\d+)\satoms    # number of atoms in group
                     """, re.VERBOSE)

def make_ndx_captured(**kwargs):
    """make_ndx that captures all output

    Standard :func:`~gromacs.make_ndx` command with the input and
    output pre-set in such a way that it can be conveniently used for
    :func:`parse_ndxlist`.

    Example::
      ndx_groups = parse_ndxlist(make_ndx_captured(n=ndx)[0])

    Note that the convenient :func:`get_ndx_groups` function does exactly
    that and can probably used in most cases.

    :Arguments:
        keywords are passed on to :func:`~gromacs.make_ndx`
    :Returns:
        (*returncode*, *output*, ``None``)
    """
    kwargs['stdout']=False   # required for proper output as described in doc
    user_input = kwargs.pop('input',[])
    user_input = [cmd for cmd in user_input if cmd != 'q']  # filter any quit
    kwargs['input'] = user_input + ['', 'q']                # necessary commands
    return gromacs.make_ndx(**kwargs)

def get_ndx_groups(ndx, **kwargs):
    """Return a list of index groups in the index file *ndx*.

    :Arguments:
        - *ndx*  is a Gromacs index file.
        - kwargs are passed to :func:`make_ndx_captured`.

    :Returns:
        list of groups as supplied by :func:`parse_ndxlist`

    Alternatively, load the index file with
    :class:`gromacs.formats.NDX` for full control.
    """
    fd, tmp_ndx = tempfile.mkstemp(suffix='.ndx')
    kwargs['o'] = tmp_ndx
    try:
        g = parse_ndxlist(make_ndx_captured(n=ndx, **kwargs)[1])
    finally:
        utilities.unlink_gmx(tmp_ndx)
    return g

def parse_ndxlist(output):
    """Parse output from make_ndx to build list of index groups::

      groups = parse_ndxlist(output)

    output should be the standard output from ``make_ndx``, e.g.::

       rc,output,junk = gromacs.make_ndx(..., input=('', 'q'), stdout=False, stderr=True)

    (or simply use

       rc,output,junk = cbook.make_ndx_captured(...)

    which presets input, stdout and stderr; of course input can be overriden.)

    :Returns:
       The function returns a list of dicts (``groups``) with fields

       name
           name of the groups
       nr
           number of the group (starts at 0)
       natoms
           number of atoms in the group

    """

    m = NDXLIST.search(output)    # make sure we pick up a proper full list
    grouplist = m.group('LIST')
    return parse_groups(grouplist)

def parse_groups(output):
    """Parse ``make_ndx`` output and return groups as a list of dicts."""
    groups = []
    for line in output.split('\n'):
        m = NDXGROUP.match(line)
        if m:
            d = m.groupdict()
            groups.append({'name': d['GROUPNAME'],
                           'nr': int(d['GROUPNUMBER']),
                           'natoms': int(d['NATOMS'])})
    return groups

class IndexBuilder(object):
    """Build an index file with specified groups and the combined group.

    This is *not* a full blown selection parser a la Charmm, VMD or
    MDAnalysis but a very quick hack.

    **Example**

       How to use the :class:`IndexBuilder`::

          G = gromacs.cbook.IndexBuilder('md_posres.pdb',
                        ['S312:OG','T313:OG1','A38:O','A309:O','@a62549 & r NA'],
                        offset=-9, out_ndx='selection.ndx')
          groupname, ndx = G.combine()
          del G

       The residue numbers are given with their canonical resids from the
       sequence or pdb. *offset=-9* says that one calculates Gromacs topology
       resids by subtracting 9 from the canonical resid.

       The combined selection is ``OR`` ed by default and written to
       *selection.ndx*. One can also add all the groups in the initial *ndx*
       file (or the :program:`make_ndx` default groups) to the output (see the
       *defaultgroups* keyword for :meth:`IndexBuilder.combine`).

       Generating an index file always requires calling
       :meth:`~IndexBuilder.combine` even if there is only a single group.

       Deleting the class removes all temporary files associated with it (see
       :attr:`IndexBuilder.indexfiles`).

    :Raises:
       If an empty group is detected (which does not always work) then a
       :exc:`gromacs.BadParameterWarning` is issued.

    :Bugs:
       If ``make_ndx`` crashes with an unexpected error then this is fairly hard to
       diagnose. For instance, in certain cases it segmentation faults when a tpr
       is provided as a *struct* file and the resulting error messages becomes ::

          GromacsError: [Errno -11] Gromacs tool failed
          Command invocation: make_ndx -o /tmp/tmp_Na1__NK7cT3.ndx -f md_posres.tpr

       In this case run the command invocation manually to see what the problem
       could be.


    .. SeeAlso:: In some cases it might be more straightforward to use
                 :class:`gromacs.formats.NDX`.

    """

    def __init__(self, struct=None, selections=None, names=None, name_all=None,
                 ndx=None, out_ndx="selection.ndx", offset=0):
        """Build a index group from the selection arguments.

        If selections and a structure file are supplied then the individual
        selections are constructed with separate calls to
        :func:`gromacs.make_ndx`. Use :meth:`IndexBuilder.combine` to combine
        them into a joint selection or :meth:`IndexBuilder.write` to simply write
        out the individual named selections (useful with *names*).

        :Arguments:

           *struct* : filename
              Structure file (tpr, pdb, ...)

           *selections* : list
              The list must contain strings or tuples, which must be be one of
              the following constructs:

                 "<1-letter aa code><resid>[:<atom name]"

                     Selects the CA of the residue or the specified atom
                     name.

                     example: ``"S312:OA"`` or ``"A22"`` (equivalent to ``"A22:CA"``)

                 ("<1-letter aa code><resid>", "<1-letter aa code><resid>, ["<atom name>"])

                    Selects a *range* of residues. If only two residue
                    identifiers are provided then all atoms are
                    selected. With an optional third atom identifier,
                    only this atom anme is selected for each residue
                    in the range. [EXPERIMENTAL]

                 "@<make_ndx selection>"

                     The ``@`` letter introduces a verbatim ``make_ndx``
                     command. It will apply the given selection without any
                     further processing or checks.

                     example: ``"@a 6234 - 6238"`` or ``'@"SOL"'`` (note the quoting)
                     or ``"@r SER & r 312 & t OA"``.

           *names* : list
              Strings to name the selections; if not supplied or if individuals
              are ``None`` then a default name is created. When simply using
              :meth:`IndexBuilder.write` then these should be supplied.

           *name_all* : string
              Name of the group that is generated by :meth:`IndexBuilder.combine`.

           *offset* : int, dict
              This number is added to the resids in the first selection scheme; this
              allows names to be the same as in a crystal structure. If offset is a
              dict then it is used to directly look up the resids.

           *ndx* : filename or list of filenames
              Optional input index file(s).

           *out_ndx* : filename
              Output index file.

        """
        self.structure = struct
        self.ndx = ndx
        self.output = out_ndx
        self.name_all = name_all
        #: *offset* as a number is added to the resids in the first selection
        #: scheme; this
        #: allows names to be the same as in a crystal structure. If *offset* is a
        #: dict then it is used to directly look up the resids. Use :meth:`gmx_resid`
        #: to transform a crystal resid to a gromacs resid.
        #:
        #: The attribute may be changed directly after init.
        self.offset = offset

        #: Auto-labelled groups use this counter.
        self._command_counter = 0

        if selections is None:
            selections = []
        if not utilities.iterable(selections):
            selections = [selections]
        self.selections = selections
        if names is None:
            names = [None] * len(selections)

        #: Specialized ``make_ndx`` that  always uses same structure
        #: and redirection (can be overridden)
        self.make_ndx = tools.Make_ndx(f=self.structure, n=self.ndx,
                                       stdout=False, stderr=False)

        #: dict, keyed by group name and pointing to index file for group
        #: (Groups are built in separate files because that is more robust
        #: as I can clear groups easily.)
        self.indexfiles = dict([self.parse_selection(selection, name)
                                for selection, name in zip(selections, names)])

    @property
    def names(self):
        """Names of all generated index groups."""
        return self.indexfiles.keys()

    def gmx_resid(self, resid):
        """Returns resid in the Gromacs index by transforming with offset."""
        try:
            gmx_resid = int(self.offset[resid])
        except (TypeError, IndexError):
            gmx_resid = resid + self.offset
        except KeyError:
            raise KeyError("offset must be a dict that contains the gmx resid for %d" % resid)
        return gmx_resid

    def combine(self, name_all=None, out_ndx=None, operation='|', defaultgroups=False):
        """Combine individual groups into a single one and write output.

        :Keywords:
           name_all : string
              Name of the combined group, ``None`` generates a name.  [``None``]
           out_ndx : filename
              Name of the output file that will contain the individual groups
              and the combined group. If ``None`` then default from the class
              constructor is used. [``None``]
           operation : character
              Logical operation that is used to generate the combined group from
              the individual groups: "|" (OR) or "&" (AND); if set to ``False``
              then no combined group is created and only the individual groups
              are written. ["|"]
           defaultgroups : bool
              ``True``: append everything to the default groups produced by
              :program:`make_ndx` (or rather, the groups provided in the ndx file on
              initialization --- if this was ``None`` then these are truly default groups);
              ``False``: only use the generated groups

        :Returns:
           ``(combinedgroup_name, output_ndx)``, a tuple showing the
           actual group name and the name of the file; useful when all names are autogenerated.

        .. Warning:: The order of the atom numbers in the combined group is
                     *not* guaranteed to be the same as the selections on input because
                     ``make_ndx`` sorts them ascending. Thus you should be careful when
                     using these index files for calculations of angles and dihedrals.
                     Use :class:`gromacs.formats.NDX` in these cases.

        .. SeeAlso:: :meth:`IndexBuilder.write`.
        """
        if not operation in ('|', '&', False):
            raise ValueError("Illegal operation %r, only '|' (OR) and '&' (AND) or False allowed." %
                             operation)
        if name_all is None:
            if operation:
                name_all = self.name_all or operation.join(self.indexfiles)
        if out_ndx is None:
            out_ndx = self.output

        if defaultgroups:
            # make a default file (using the original ndx where provided!!)
            fd, default_ndx = tempfile.mkstemp(suffix='.ndx', prefix='default__')
            try:
                self.make_ndx(o=default_ndx, input=['q'])
            except:
                utilities.unlink_gmx(default_ndx)
                raise
            ndxfiles = [default_ndx]
        else:
            ndxfiles = []

        ndxfiles.extend(self.indexfiles.values())

        if operation:
            # combine multiple selections and name them
            try:
                fd, tmp_ndx = tempfile.mkstemp(suffix='.ndx', prefix='combined__')
                # combine all selections by loading ALL temporary index files
                operation = ' '+operation.strip()+' '
                cmd = [operation.join(['"%s"' % gname for gname in self.indexfiles]),
                       '', 'q']
                rc,out,err = self.make_ndx(n=ndxfiles, o=tmp_ndx, input=cmd)
                if self._is_empty_group(out):
                    warnings.warn("No atoms found for %(cmd)r" % vars(),
                                  category=BadParameterWarning)

                # second pass for naming, sigh (or: use NDX ?)
                groups = parse_ndxlist(out)
                last = groups[-1]
                # name this group
                name_cmd = ["name %d %s" % (last['nr'], name_all), 'q']
                rc,out,err = self.make_ndx(n=tmp_ndx, o=out_ndx, input=name_cmd)
                # For debugging, look at out and err or set stdout=True, stderr=True
                # TODO: check out if at least 1 atom selected
                ##print "DEBUG: combine()"
                ##print out
            finally:
                utilities.unlink_gmx(tmp_ndx)
                if defaultgroups:
                    utilities.unlink_gmx(default_ndx)
        else:
            # just write individual groups in one file (name_all --> None)
            rc,out,err = self.make_ndx(n=ndxfiles, o=out_ndx, input=['','q'])

        return name_all, out_ndx

    def write(self, out_ndx=None, defaultgroups=False):
        """Write individual (named) groups to *out_ndx*."""
        name_all, out_ndx = self.combine(operation=False, out_ndx=out_ndx, defaultgroups=defaultgroups)
        return out_ndx

    def cat(self, out_ndx=None):
        """Concatenate input index files.

        Generate a new index file that contains the default Gromacs index
        groups (if a structure file was defined) and all index groups from the
        input index files.

        :Arguments:
           out_ndx : filename
              Name of the output index file; if ``None`` then use the default
              provided to the constructore. [``None``].
        """
        if out_ndx is None:
            out_ndx = self.output
        self.make_ndx(o=out_ndx, input=['q'])
        return out_ndx

    def parse_selection(self, selection, name=None):
        """Retuns (groupname, filename) with index group."""

        if type(selection) is tuple:
            # range
            process = self._process_range
        elif selection.startswith('@'):
            # verbatim make_ndx command
            process = self._process_command
            selection = selection[1:]
        else:
            process = self._process_residue
        return process(selection, name)

    def _process_command(self, command, name=None):
        """Process ``make_ndx`` command and  return name and temp index file."""

        self._command_counter += 1
        if name is None:
            name = "CMD%03d" % self._command_counter

        # Need to build it with two make_ndx calls because I cannot reliably
        # name the new group without knowing its number.
        try:
            fd, tmp_ndx = tempfile.mkstemp(suffix='.ndx', prefix='tmp_'+name+'__')
            cmd = [command, '', 'q']   # empty command '' necessary to get list
            # This sometimes fails with 'OSError: Broken Pipe' --- hard to debug
            rc,out,err = self.make_ndx(o=tmp_ndx, input=cmd)
            self.check_output(out, "No atoms found for selection %(command)r." % vars(), err=err)
            # For debugging, look at out and err or set stdout=True, stderr=True
            # TODO: check '  0 r_300_&_ALA_&_O     :     1 atoms' has at least 1 atom
            ##print "DEBUG: _process_command()"
            ##print out
            groups = parse_ndxlist(out)
            last = groups[-1]
            # reduce and name this group
            fd, ndx = tempfile.mkstemp(suffix='.ndx', prefix=name+'__')
            name_cmd = ["keep %d" % last['nr'],
                        "name 0 %s" % name, 'q']
            rc,out,err = self.make_ndx(n=tmp_ndx, o=ndx, input=name_cmd)
        finally:
            utilities.unlink_gmx(tmp_ndx)

        return name, ndx

    #: regular expression to match and parse a residue-atom selection
    RESIDUE = re.compile("""
                 (?P<aa>([ACDEFGHIKLMNPQRSTVWY])   # 1-letter amino acid
                        |                          #   or
                        ([A-Z][A-Z][A-Z][A-Z]?)    # 3-letter or 4-letter residue name
                 )
                 (?P<resid>\d+)                    # resid
                 (:                                # separator ':'
                   (?P<atom>\w+)                   # atom name
                 )?                                # possibly one
            """, re.VERBOSE | re.IGNORECASE)

    def _process_residue(self, selection, name=None):
        """Process residue/atom selection and return name and temp index file."""

        if name is None:
            name = selection.replace(':', '_')

        # XXX: use _translate_residue() ....
        m = self.RESIDUE.match(selection)
        if not m:
            raise ValueError("Selection %(selection)r is not valid." % vars())

        gmx_resid = self.gmx_resid(int(m.group('resid')))
        residue = m.group('aa')
        if len(residue) == 1:
            gmx_resname = utilities.convert_aa_code(residue) # only works for AA
        else:
            gmx_resname = residue                            # use 3-letter for any resname
        gmx_atomname = m.group('atom')
        if gmx_atomname is None:
            gmx_atomname = 'CA'

        #: select residue <gmx_resname><gmx_resid> atom <gmx_atomname>
        _selection = 'r %(gmx_resid)d & r %(gmx_resname)s & a %(gmx_atomname)s' % vars()
        cmd = ['keep 0', 'del 0',
               _selection,
               'name 0 %(name)s' % vars(),
               'q']
        fd, ndx = tempfile.mkstemp(suffix='.ndx', prefix=name+'__')
        rc,out,err = self.make_ndx(n=self.ndx, o=ndx, input=cmd)
        self.check_output(out, "No atoms found for "
                          "%(selection)r --> %(_selection)r" % vars())
        # For debugging, look at out and err or set stdout=True, stderr=True
        ##print "DEBUG: _process_residue()"
        ##print out

        return name, ndx

    def _process_range(self, selection, name=None):
        """Process a range selection.

        ("S234", "A300", "CA")   --> selected all CA in this range
        ("S234", "A300")         --> selected all atoms in this range

        .. Note:: Ignores residue type, only cares about the resid (but still required)
        """

        try:
            first, last, gmx_atomname = selection
        except ValueError:
            try:
                first, last = selection
                gmx_atomname = '*'
            except:
                logger.error("%r is not a valid range selection", selection)
                raise
        if name is None:
            name = "%(first)s-%(last)s_%(gmx_atomname)s" % vars()

        _first = self._translate_residue(first, default_atomname=gmx_atomname)
        _last = self._translate_residue(last, default_atomname=gmx_atomname)

        _selection = 'r %d - %d & & a %s' % (_first['resid'],  _last['resid'], gmx_atomname)
        cmd = ['keep 0', 'del 0',
               _selection,
               'name 0 %(name)s' % vars(),
               'q']
        fd, ndx = tempfile.mkstemp(suffix='.ndx', prefix=name+'__')
        rc,out,err = self.make_ndx(n=self.ndx, o=ndx, input=cmd)
        self.check_output(out, "No atoms found for "
                          "%(selection)r --> %(_selection)r" % vars())
        # For debugging, look at out and err or set stdout=True, stderr=True
        ##print "DEBUG: _process_residue()"
        ##print out

        return name, ndx


    def _translate_residue(self, selection, default_atomname='CA'):
        """Translate selection for a single res to make_ndx syntax."""
        m = self.RESIDUE.match(selection)
        if not m:
            errmsg = "Selection %(selection)r is not valid." % vars()
            logger.error(errmsg)
            raise ValueError(errmsg)

        gmx_resid = self.gmx_resid(int(m.group('resid')))    # magic offset correction
        residue = m.group('aa')
        if len(residue) == 1:
            gmx_resname = utilities.convert_aa_code(residue) # only works for AA
        else:
            gmx_resname = residue                            # use 3-letter for any resname

        gmx_atomname = m.group('atom')
        if gmx_atomname is None:
            gmx_atomname = default_atomname

        return {'resname':gmx_resname, 'resid':gmx_resid, 'atomname':gmx_atomname}



    def check_output(self, make_ndx_output, message=None, err=None):
        """Simple tests to flag problems with a ``make_ndx`` run."""
        if message is None:
            message = ""
        else:
            message = '\n' + message
        def format(output, w=60):
            hrule = "====[ GromacsError (diagnostic output) ]".ljust(w,"=")
            return hrule + '\n' + str(output) + hrule

        rc = True
        if self._is_empty_group(make_ndx_output):
            warnings.warn("Selection produced empty group.%(message)s"
                          % vars(), category=GromacsValueWarning)
            rc = False
        if self._has_syntax_error(make_ndx_output):
            rc = False
            out_formatted = format(make_ndx_output)
            raise GromacsError("make_ndx encountered a Syntax Error, "
                               "%(message)s\noutput:\n%(out_formatted)s" % vars())
        if make_ndx_output.strip() == "":
            rc = False
            out_formatted = format(err)
            raise GromacsError("make_ndx produced no output, "
                               "%(message)s\nerror output:\n%(out_formatted)s" % vars())
        return rc

    def _is_empty_group(self, make_ndx_output):
        m = re.search('Group is empty', make_ndx_output)
        return not (m is None)

    def _has_syntax_error(self, make_ndx_output):
        m = re.search('Syntax error:', make_ndx_output)
        return not (m is None)

    def __del__(self):
        try:
            for path in self.indexfiles.values():
                utilities.unlink_gmx(path)
                # Removes auto-backup files, too (which we have because mkstemp creates
                # an empty file and make_ndx backs that up).
        except (AttributeError, OSError):
            # all exceptions are ignored inside __del__ anyway but these
            # two we do not even want to be noticed off:
            # AttributeError: when reloading the module, OSError: when file disappeared
            pass


class Transformer(utilities.FileUtils):
    """Class to handle transformations of trajectories.

    1. Center, compact, and fit to reference structure in tpr
       (optionally, only center in the xy plane): :meth:`~Transformer.center_fit`
    2. Write compact xtc and tpr with water removed: :meth:`~Transformer.strip_water`
    3. Write compact xtc and tpr only with protein: :meth:`~Transformer.keep_protein_only`

    """

    def __init__(self, s="topol.tpr", f="traj.xtc", n=None, force=None, 
                 dirname=os.path.curdir, outdir=None):
        """Set up Transformer with structure and trajectory.

        Supply *n* = tpr, *f* = xtc (and *n* = ndx) relative to dirname.

        :Keywords:
           *s*
              tpr file (or similar); note that this should not contain
              position restraints if it is to be used with a reduced
              system (see :meth:`~Transformer.strip_water`)
           *f*
              trajectory (xtc, trr, ...)
           *n*
              index file (it is typically safe to leave this as ``None``; in
              cases where a trajectory needs to be centered on non-standard
              groups this should contain those groups)
           *force*
              Set the default behaviour for handling existing files:
                - ``True``: overwrite existing trajectories
                - ``False``: throw a IOError exception
                - ``None``: skip existing and log a warning [default]
           *dirname*
              directory in which all operations are performed, relative paths
              are interpreted relative to *dirname* [.]
           *outdir*
              directory under which output files are placed; by default
              the same directory where the input files live
        """

        self.tpr = self.filename(s, ext="tpr", use_my_ext=True)
        self.xtc = self.filename(f, ext="xtc", use_my_ext=True)
        self.ndx = n
        self.dirname = dirname
        self.outdir = utilities.realpath(outdir) if outdir is not None else None
        self.force = force
        self.nowater = {}     # data for trajectory stripped from water
        self.proteinonly = {} # data for a protein-only trajectory

        with utilities.in_dir(self.dirname, create=False):
            for f in (self.tpr, self.xtc, self.ndx):
                if f is None:
                    continue
                if not os.path.exists(f):
                    msg = "Possible problem: File %(f)r not found in %(dirname)r." % vars()
                    warnings.warn(msg, category=MissingDataWarning)
                    logger.warn(msg)
        logger.info("%r initialised", self)

    def __repr__(self):
            return "%s(s=%r, f=%r, n=%r, force=%r)" % (self.__class__.__name__,
                                                       self.tpr, self.xtc, self.ndx, self.force)
    def outfile(self, p):
        """Path for an output file.

        If :attr:`outdir` is set then the path is
        ``outdir/basename(p)`` else just ``p``
        """
        if self.outdir is not None:
            return os.path.join(self.outdir, os.path.basename(p))
        else:
            return p

    def rp(self, *args):
        """Return canonical path to file under *dirname* with components *args*

         If *args* form an absolute path then just return it as the absolute path.
         """
        try:
            p = os.path.join(*args)
            if os.path.isabs(p):
                return p
        except TypeError:
            pass
        return utilities.realpath(self.dirname, *args)

    def center_fit(self, **kwargs):
        """Write compact xtc that is fitted to the tpr reference structure.

        See :func:`gromacs.cbook.trj_fitandcenter` for details and
        description of *kwargs* (including *input*, *input1*, *n* and 
        *n1* for how to supply custom index groups). The most important ones are listed
        here but in most cases the defaults should work.

        :Keywords:
           *s*
             Input structure (typically the default tpr file but can be set to
             some other file with a different conformation for fitting)
           *n*
             Alternative index file.
           *o*
             Name of the output trajectory.
           *xy* : Boolean
             If ``True`` then only fit in xy-plane (useful for a membrane normal
             to z). The default is ``False``.
           *force*
             - ``True``: overwrite existing trajectories
             - ``False``: throw a IOError exception
             - ``None``: skip existing and log a warning [default]

        :Returns:
              dictionary with keys *tpr*, *xtc*, which are the names of the
              the new files
        """
        kwargs.setdefault('s', self.tpr)
        kwargs.setdefault('n', self.ndx)
        kwargs['f'] = self.xtc
        kwargs.setdefault('o', self.outfile(self.infix_filename(None, self.xtc, '_centfit', 'xtc')))
        force = kwargs.pop('force', self.force)

        logger.info("Centering and fitting trajectory %(f)r..." % kwargs)
        with utilities.in_dir(self.dirname):
            if not self.check_file_exists(kwargs['o'], resolve="indicate", force=force):
                trj_fitandcenter(**kwargs)
            logger.info("Centered and fit trajectory: %(o)r." % kwargs)
        return {'tpr': self.rp(kwargs['s']), 'xtc': self.rp(kwargs['o'])}

    def fit(self, xy=False, **kwargs):
        """Write xtc that is fitted to the tpr reference structure.

        Runs :class:`gromacs.tools.trjconv` with appropriate arguments
        for fitting. The most important *kwargs* are listed
        here but in most cases the defaults should work.

        Note that the default settings do *not* include centering or
        periodic boundary treatment as this often does not work well
        with fitting. It is better to do this as a separate step (see
        :meth:`center_fit` or :func:`gromacs.cbook.trj_fitandcenter`)

        :Keywords:
           *s*
             Input structure (typically the default tpr file but can be set to
             some other file with a different conformation for fitting)
           *n*
             Alternative index file.
           *o*
             Name of the output trajectory. A default name is created.
             If e.g. *dt* = 100  is one of the *kwargs* then the default name includes
             "_dt100ps".
          *xy* : boolean
             If ``True`` then only do a rot+trans fit in the xy plane
             (good for membrane simulations); default is ``False``.
          *force*
            ``True``: overwrite existing trajectories
            ``False``: throw a IOError exception
            ``None``: skip existing and log a warning [default]
          *fitgroup*
            index group to fit on ["backbone"]

            .. Note:: If keyword *input* is supplied then it will override
                      *fitgroup*; *input* = ``[fitgroup, outgroup]``
          *kwargs*
             kwargs are passed to :func:`~gromacs.cbook.trj_xyfitted`

        :Returns:
              dictionary with keys *tpr*, *xtc*, which are the names of the
              the new files
        """
        kwargs.setdefault('s', self.tpr)
        kwargs.setdefault('n', self.ndx)
        kwargs['f'] = self.xtc
        force = kwargs.pop('force', self.force)
        if xy:
            fitmode = 'rotxy+transxy'
            kwargs.pop('fit', None)
            infix_default = '_fitxy'
        else:
            fitmode = kwargs.pop('fit', 'rot+trans')  # user can use 'progressive', too
            infix_default = '_fit'

        dt = kwargs.get('dt', None)
        if dt:
            infix_default += '_dt%dps' % int(dt)    # dt in ps

        kwargs.setdefault('o', self.outfile(self.infix_filename(None, self.xtc, infix_default, 'xtc')))
        fitgroup = kwargs.pop('fitgroup', 'backbone')
        kwargs.setdefault('input', [fitgroup, "system"])

        if kwargs.get('center', False):
            logger.warn("Transformer.fit(): center=%(center)r used: centering should not be combined with fitting.", kwargs)
            if len(kwargs['inputs']) != 3:
                logger.error("If you insist on centering you must provide three groups in the 'input' kwarg: (center, fit, output)")
                raise ValuError("Insufficient index groups for centering,fitting,output")

        logger.info("Fitting trajectory %r to with xy=%r...", kwargs['f'], xy)
        logger.info("Fitting on index group %(fitgroup)r", vars())
        with utilities.in_dir(self.dirname):
            if self.check_file_exists(kwargs['o'], resolve="indicate", force=force):
                logger.warn("File %r exists; force regenerating it with force=True.", kwargs['o'])
            else:
                gromacs.trjconv(fit=fitmode, **kwargs)
                logger.info("Fitted trajectory (fitmode=%s): %r.", fitmode, kwargs['o'])
        return {'tpr': self.rp(kwargs['s']), 'xtc': self.rp(kwargs['o'])}

    def strip_water(self, os=None, o=None, on=None, compact=False, 
                    resn="SOL", groupname="notwater", **kwargs):
        """Write xtc and tpr with water (by resname) removed.

        :Keywords:
           *os*
              Name of the output tpr file; by default use the original but
              insert "nowater" before suffix.
           *o*
              Name of the output trajectory; by default use the original name but
              insert "nowater" before suffix.
           *on*
              Name of a new index file (without water).
           *compact*
              ``True``: write a compact and centered trajectory
              ``False``: use trajectory as it is [``False``]
           *centergroup*
              Index group used for centering ["Protein"]

              .. Note:: If *input* is provided (see below under *kwargs*)
                        then *centergroup* is ignored and the group for 
                        centering is taken as the first entry in *input*.

           *resn*
              Residue name of the water molecules; all these residues are excluded.
           *groupname*
              Name of the group that is generated by subtracting all waters
              from the system.
           *force* : Boolean
             - ``True``: overwrite existing trajectories
             - ``False``: throw a IOError exception
             - ``None``: skip existing and log a warning [default]
           *kwargs*
              are passed on to :func:`gromacs.cbook.trj_compact` (unless the
              values have to be set to certain values such as s, f, n, o
              keywords). The *input* keyword is always mangled: Only the first
              entry (the group to centre the trajectory on) is kept, and as a
              second group (the output group) *groupname* is used.

        :Returns:
              dictionary with keys *tpr*, *xtc*, *ndx* which are the names of the
              the new files

        .. warning:: The input tpr file should *not* have *any position restraints*;
                     otherwise Gromacs will throw a hissy-fit and say

                     *Software inconsistency error: Position restraint coordinates are
                     missing*

                     (This appears to be a bug in Gromacs 4.x.)
        """
        force = kwargs.pop('force', self.force)

        newtpr = self.outfile(self.infix_filename(os, self.tpr, '_nowater'))
        newxtc = self.outfile(self.infix_filename(o, self.xtc, '_nowater'))
        newndx = self.outfile(self.infix_filename(on, self.tpr, '_nowater', 'ndx'))

        nowater_ndx = self._join_dirname(newtpr, "nowater.ndx")    # refers to original tpr

        if compact:
            TRJCONV = trj_compact
            # input overrides centergroup
            if kwargs.get('centergroup') is not None and 'input' in kwargs:
		logger.warn("centergroup = %r will be superceded by input[0] = %r", kwargs['centergroup'], kwargs['input'][0])
            _input = kwargs.get('input', [kwargs.get('centergroup', 'Protein')])
            kwargs['input'] = [_input[0], groupname]  # [center group, write-out selection]
            del _input
	    logger.info("Creating a compact trajectory centered on group %r", kwargs['input'][0])
	    logger.info("Writing %r to the output trajectory", kwargs['input'][1])
        else:
            TRJCONV = gromacs.trjconv
            kwargs['input'] = [groupname]
            logger.info("Writing %r to the output trajectory (no centering)", kwargs['input'][0])
        # clean kwargs, only legal arguments for Gromacs tool trjconv should remain
        kwargs.pop("centergroup", None)

        NOTwater = "! r %(resn)s" % vars()  # make_ndx selection ("not water residues")
        with utilities.in_dir(self.dirname):
            # ugly because I cannot break from the block
            if not self.check_file_exists(newxtc, resolve="indicate", force=force):
                # make no-water index
                B = IndexBuilder(struct=self.tpr, selections=['@'+NOTwater],
                                 ndx=self.ndx, out_ndx=nowater_ndx)
                B.combine(name_all=groupname, operation="|", defaultgroups=True)
                logger.debug("Index file for water removal: %r", nowater_ndx)

                logger.info("TPR file without water %(newtpr)r" % vars())
                gromacs.tpbconv(s=self.tpr, o=newtpr, n=nowater_ndx, input=[groupname])

                logger.info("NDX of the new system %r", newndx)
                gromacs.make_ndx(f=newtpr, o=newndx, input=['q'], stderr=False, stdout=False)
		# PROBLEM: If self.ndx contained a custom group required for fitting then we are loosing
                #          this group here. We could try to merge only this group but it is possible that
                #          atom indices changed. The only way to solve this is to regenerate the group with
                #          a selection or only use Gromacs default groups.

                logger.info("Trajectory without water %(newxtc)r" % vars())
                kwargs['s'] = self.tpr
                kwargs['f'] = self.xtc
                kwargs['n'] = nowater_ndx
                kwargs['o'] = newxtc
                TRJCONV(**kwargs)

                logger.info("pdb and gro for visualization")
                for ext in 'pdb', 'gro':
                    try:
                        # see warning in doc ... so we don't use the new xtc but the old one
                        kwargs['o'] = self.filename(newtpr, ext=ext)
                        TRJCONV(dump=0, stdout=False, stderr=False, **kwargs)  # silent
                    except:
                        logger.exception("Failed building the water-less %(ext)s. "
                                         "Position restraints in tpr file (see docs)?" % vars())
            logger.info("strip_water() complete")

        self.nowater[self.rp(newxtc)] = Transformer(dirname=self.dirname, s=newtpr,
                                           f=newxtc, n=newndx, force=force)
        return {'tpr':self.rp(newtpr), 'xtc':self.rp(newxtc), 'ndx':self.rp(newndx)}


    # TODO: could probably unify strip_water() and keep_protein_only()
    # (given that the latter was produced by copy&paste+search&replace...)

    def keep_protein_only(self, os=None, o=None, on=None, compact=False,
                          groupname="proteinonly", **kwargs):
        """Write xtc and tpr only containing the protein.

        :Keywords:
           *os*
              Name of the output tpr file; by default use the original but
              insert "proteinonly" before suffix.
           *o*
              Name of the output trajectory; by default use the original name but
              insert "proteinonly" before suffix.
           *on*
              Name of a new index file.
           *compact*
              ``True``: write a compact and centered trajectory
              ``False``: use trajectory as it is [``False``]
           *groupname*
              Name of the protein-only group.
           *keepalso*
              List of literal make_ndx selections of additional groups that should
              be kept, e.g. ['resname DRUG', 'atom 6789'].
           *force* : Boolean
             - ``True``: overwrite existing trajectories
             - ``False``: throw a IOError exception
             - ``None``: skip existing and log a warning [default]
           *kwargs*
              are passed on to :func:`gromacs.cbook.trj_compact` (unless the
              values have to be set to certain values such as s, f, n, o
              keywords). The *input* keyword is always mangled: Only the first
              entry (the group to centre the trajectory on) is kept, and as a
              second group (the output group) *groupname* is used.

        :Returns:
              dictionary with keys *tpr*, *xtc*, *ndx* which are the names of the
              the new files

        .. warning:: The input tpr file should *not* have *any position restraints*;
                     otherwise Gromacs will throw a hissy-fit and say

                     *Software inconsistency error: Position restraint coordinates are
                     missing*

                     (This appears to be a bug in Gromacs 4.x.)
        """
        force = kwargs.pop('force', self.force)
        suffix = 'proteinonly'
        newtpr = self.outfile(self.infix_filename(os, self.tpr, '_'+suffix))
        newxtc = self.outfile(self.infix_filename(o, self.xtc, '_'+suffix))
        newndx = self.outfile(self.infix_filename(on, self.tpr, '_'+suffix, 'ndx'))

        selection_ndx = suffix+".ndx"    # refers to original tpr

        if compact:
            TRJCONV = trj_compact
            _input = kwargs.get('input', ['Protein'])
            kwargs['input'] = [_input[0], groupname]  # [center group, write-out selection]
            del _input
        else:
            TRJCONV = gromacs.trjconv
            kwargs['input'] = [groupname]

        selections = ['@'+sel for sel in ['"Protein"'] + kwargs.pop('keepalso',[])]
        with utilities.in_dir(self.dirname):
            # ugly because I cannot break from the block
            if not self.check_file_exists(newxtc, resolve="indicate", force=force):
                # make index (overkill for 'Protein' but maybe we want to enhance
                # it in the future, e.g. with keeping ions/ligands as well?
                B = IndexBuilder(struct=self.tpr, selections=selections,
                                 ndx=self.ndx, out_ndx=selection_ndx)
                B.combine(name_all=groupname, operation="|", defaultgroups=True)

                logger.info("TPR file containg the protein %(newtpr)r" % vars())
                gromacs.tpbconv(s=self.tpr, o=newtpr, n=selection_ndx, input=[groupname])

                logger.info("NDX of the new system %(newndx)r" % vars())
                gromacs.make_ndx(f=newtpr, o=newndx, input=['q'], stderr=False, stdout=False)

                logger.info("Trajectory with only the protein %(newxtc)r" % vars())
                kwargs['s'] = self.tpr
                kwargs['f'] = self.xtc
                kwargs['n'] = selection_ndx
                kwargs['o'] = newxtc
                TRJCONV(**kwargs)

                logger.info("pdb and gro for visualization")
                for ext in 'pdb', 'gro':
                    try:
                        # see warning in doc ... so we don't use the new xtc but the old one
                        kwargs['o'] = self.filename(newtpr, ext=ext)
                        TRJCONV(dump=0, stdout=False, stderr=False, **kwargs)  # silent
                    except:
                        logger.exception("Failed building the protein-only %(ext)s. "
                                         "Position restraints in tpr file (see docs)?" % vars())
            logger.info("keep_protein_only() complete")

        self.proteinonly[self.rp(newxtc)] = Transformer(dirname=self.dirname, s=newtpr,
                                               f=newxtc, n=newndx, force=force)
        return {'tpr':self.rp(newtpr), 'xtc':self.rp(newxtc), 'ndx':self.rp(newndx)}

    def strip_fit(self, **kwargs):
        """Strip water and fit to the remaining system.

        First runs :meth:`strip_water` and then :meth:`fit`; see there
        for arguments.

        - *strip_input* is used for :meth:`strip_water` (but is only useful in
          special cases, e.g. when there is no Protein group defined. Then set
          *strip_input* = ``['Other']``.

        - *input* is passed on to :meth:`fit` and can contain the
          ``[center_group, fit_group, output_group]``

        - *fitgroup* is only passed to :meth:`fit` and just contains
          the group to fit to ("backbone" by default)

	  .. warning:: *fitgroup* can only be a Gromacs default group and not
                       a custom group (because the indices change after stripping)

        - By default *fit* = "rot+trans" (and *fit* is passed to :meth:`fit`,
          together with the *xy* = ``False`` keyword)

        .. Note:: The call signature of :meth:`strip_water` is somewhat different from this one.
        """
        kwargs.setdefault('fit', 'rot+trans')
        kw_fit = {}
        for k in ('xy', 'fit', 'fitgroup', 'input'):
            if k in kwargs:
                kw_fit[k] = kwargs.pop(k)

        kwargs['input'] = kwargs.pop('strip_input', ['Protein'])
        kwargs['force'] = kw_fit['force'] = kwargs.pop('force', self.force)

        paths = self.strip_water(**kwargs)    # updates self.nowater
        transformer_nowater = self.nowater[paths['xtc']]  # make sure to get the one we just produced
        return transformer_nowater.fit(**kw_fit)          # use new Transformer's fit()

    def _join_dirname(self, *args):
        """return os.path.join(os.path.dirname(args[0]), *args[1:])"""
        # extra function because I need to use it in a method that defines
        # the kwarg 'os', which collides with os.path...
        return os.path.join(os.path.dirname(args[0]), *args[1:])



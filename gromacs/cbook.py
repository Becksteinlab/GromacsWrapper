# $Id$
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
problem (see also `Manipulating trajectories`_):

.. function:: rmsd_backbone([s="md.tpr", f="md.xtc"[, ...]])

   Computes the RMSD of the "Backbone" atoms after fitting to the
   "Backbone" (including both translation and rotation).


Manipulating trajectories
-------------------------

Standard invocations for compacting or fitting trajectories.

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
.. autoclass:: Transformer
   :members:


Processing output
-----------------

There are cases when a script has to to do different things depending
on the output from a Gromacs tool. 

For instance, a common case is to check the total charge after
grompping a tpr file. The ``grompp_qtot`` function does just that.

.. autofunction:: grompp_qtot
.. autofunction:: parse_ndxlist


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

from __future__ import with_statement

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
from gromacs import GromacsError, BadParameterWarning
import tools
import utilities

trj_compact = tools.Trjconv(ur='compact', center=True, boxcenter='tric', pbc='mol',
                            input=('protein','system'),
                            doc="""
Writes a compact representation of the system centered on the protein""")

rmsd_backbone = tools.G_rms(what='rmsd', fit='rot+trans',
                            input=('Backbone','Backbone'),
                            doc="""
Computes RMSD of backbone after fitting to the backbone.""")

# Gromacs 4.x
trj_xyfitted = tools.Trjconv(fit='rotxy+transxy',
                            center=True, boxcenter='rect', pbc='whole',
                            input=('backbone', 'protein','system'),
                            doc="""
Writes a trajectory centered and fitted to the protein in the XY-plane only.

This is useful for membrane proteins. The system *must* be oriented so
that the membrane is in the XY plane. The protein backbone is used
for the least square fit, centering is done for the whole protein.

.. Note:: Gromacs 4.x only""")

def trj_fitandcenter(xy=False, **kwargs):
    """Center everything and make a compact representation (pass 1) and fit the system to a reference (pass 2).

    :Keywords:
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
    *dump*, *timestep*; in this case do things manually.

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

       trjconv -s a.tpr -f a.xtc -o b.xtc -center tric -ur compact -pbc none
       trjconv -s a.tpr -f b.xtc -o c.xtc -fit rot+trans
    
    .. _`g_spatial documentation`: http://oldwiki.gromacs.org/index.php/Manual:g_spatial_4.0.3
    """
    if xy:
        fitmode = 'rotxy+transxy'
    else:
        fitmode = 'rot+trans'
        
    intrj = kwargs.pop('f', None)
    # get the correct suffix for the intermediate step: only trr will
    # keep velocities/forces!
    suffix = os.path.splitext(intrj)[1]
    if not suffix in ('xtc', 'trr'):
        suffix = '.xtc'
    outtrj = kwargs.pop('o', None)    
    inpfit = kwargs.pop('input', ('backbone', 'protein','system'))
    try:
        inpcompact = inpfit[1:]     # use 2nd and 3rd group for compact
    except TypeError:
        inpcompact = None
    fd, tmptrj = tempfile.mkstemp(suffix=suffix, prefix='fitted_')

    logger.info("Input trajectory:  %(intrj)r\nOutput trajectory: %(outtrj)r"% vars())
    logger.info("... writing temporary trajectory %(tmptrj)r (will be auto-cleaned)." % vars())
    sys.stdout.flush()
    try:
        trj_compact(f=intrj, o=tmptrj, input=inpcompact, **kwargs)        
        trj_xyfitted(f=tmptrj, o=outtrj, fit=fitmode, input=inpfit, **kwargs)
    finally:
        utilities.unlink_gmx(tmptrj)


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

def grompp_qtot(*args, **kwargs):
    """Run ``gromacs.grompp`` and return the total charge of the  system.

    :Arguments:  
       The arguments are the ones one would pass to :func:`gromacs.grompp`.
    :Returns:
       The total charge as reported

    .. note::

       * The stdout output of grompp is not shown. This can make debugging
         pretty hard.  Try running the normal :func:`gromacs.grompp` command and
         analyze the output if the debugging messages are not sufficient.
       * Check that ``qtot`` is correct; because the function is based on pattern 
         matching of the output it can break when the output format changes.

    """

    # match '  System has non-zero total charge: -4.000001e+00',
    qtot_pattern = re.compile('System has non-zero total charge: *(?P<qtot>[-+]?\d*\.\d+[eE][-+]\d+)')
    kwargs['stdout'] = False
    rc, output, junk = gromacs.grompp(*args, **kwargs)
    qtot = 0
    for line in output.split('\n'):
        m = qtot_pattern.search(line)
        if m:
            qtot = float(m.group('qtot'))
            break
    logger.info("system total charge qtot = %(qtot)r" % vars())
    return qtot


# Editing textual input files
# ---------------------------

def edit_mdp(mdp, new_mdp=None, **substitutions):
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
        *substitutions*
            parameter=value pairs, where parameter is defined by the Gromacs mdp file; 
            dashes in parameter names have to be replaced by underscores.

    :Returns:    
        Dict of parameters that have *not* been substituted.

    Example::

       edit_mdp('md.mdp', new_mdp='long_md.mdp', nsteps=100000, nstxtcout=1000, lincs_iter=2)

    .. Note::
    
       * Dashes in Gromacs mdp parameters have to be replaced by an underscore
         when supplied as python keyword arguments (a limitation of python). For example
         the MDP syntax is  ``lincs-iter = 4`` but the corresponding  keyword would be 
         ``lincs_iter = 4``.
       * If the keyword is set as a dict key, eg ``mdp_params['lincs-iter']=4`` then one
         does not have to substitute.
       * Parameters *aa_bb* and *aa-bb* are considered the same (although this should not be 
         a problem in practice because there are no mdp parameters that only differ by a underscore).
       * This code is more compact in ``Perl`` as one can use ``s///`` operators:
         ``s/^(\s*${key}\s*=\s*).*/$1${val}/``

    """
    if new_mdp is None:
        new_mdp = mdp

    # None parameters should be ignored (simple way to keep the template defaults)
    substitutions = dict([(k,v) for k,v in substitutions.items() if not v is None])

    params = substitutions.keys()[:]   # list will be reduced for each match

    def demangled(p):
        """Return a RE string that matches the parameter."""
        return p.replace('_', '[-_]')  # must catch either - or _

    patterns = dict([(parameter,
                      re.compile("""\
                       (?P<assignment>\s*%s\s*=\s*)  # parameter == everything before the value
                       (?P<value>[^;]*)              # value (stop before comment=;)
                       (?P<comment>\s*;.*)?          # optional comment           
                       """ % demangled(parameter), re.VERBOSE))
                     for parameter in substitutions])

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
                    new_line = assignment + str(substitutions[p]) + comment
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
    return dict([(p, substitutions[p]) for p in params])

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

    .. note::
    
       * No sanity checks are performed and the substitutions must be supplied
         exactly as shown.    
       * Only the first matching substitution is applied to a line; thus the order
         of the substitution commands matters. This behaviour was chosen to avoid
         ambiguity when a substitution would create a match for a subsequent rule.
    """
    if newname is None:
        newname = filename

    # no sanity checks (figure out later how to give decent diagnostics)
    _substitutions = [{'lRE': re.compile(str(lRE)),
                       'sRE': re.compile(str(sRE)),
                       'repl': str(repl)}
                      for lRE,sRE,repl in substitutions]

    target = tempfile.TemporaryFile()
    with open(filename) as src:
        logger.info("editing txt = %r (%d substitutions)" % (filename, len(substitutions)))
        for line in src:
            new_line = line[:]
            for subst in _substitutions:
                m = subst['lRE'].match(line)    
                if m:              # apply substition to this line?
                    logger.debug('match:    '+line)
                    new_line = subst['sRE'].sub(subst['repl'], line)
                    logger.debug('replaced: '+new_line)
                    break   # only apply the first matching substitution!
            target.write(new_line)
    target.seek(0)

    with open(newname, 'w') as final:
        shutil.copyfileobj(target, final)
    target.close()
    logger.info("edited txt = %(newname)r" % vars())


# Working with index files and index groups
# -----------------------------------------

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
#: compiler regular expression to match a single line of 
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
    kwargs['stderr']=True    # ...
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
        them into a joint selection.

        :Arguments:

           struct : filename
              Structure file (tpr, pdb, ...)

           selections : list
              The list must contain strings, which must be be one of
              the following constructs:

                 "<1-letter aa code><resid>[:<atom name]"

                     Selects the CA of the residue or the specified atom
                     name.

                     example: ``"S312:OA"`` or ``"A22"`` (equivalent to ``"A22:CA"``)

                 "@<make_ndx selection>"

                     The ``@`` letter introduces a verbatim ``make_ndx``
                     command. It will apply the given selection without any
                     further processing or checks.

                     example: ``"@a 6234 - 6238"`` or ``'@"SOL"'`` (note the quoting)
                     or ``"@r SER & r 312 & t OA"``.

           names : list
              Strings to name the selections; if not supplied or if individuals
              are ``None`` then a default name is created.

           offset : int, dict
              This number is added to the resids in the first selection scheme; this
              allows names to be the same as in a crystal structure. If offset is a 
              dict then it is used to directly look up the resids.

           ndx : filename or list of filenames
              Optional input index file(s).

           out_ndx : filename
              Output index file.
              
        """
        self.structure = struct
        self.ndx = ndx
        self.output = out_ndx
        self.name_all = name_all
        #: This number is added to the resids in the first selection scheme; this
        #: allows names to be the same as in a crystal structure. If offset is a 
        #: dict then it is used to directly look up the resids. Use :meth:`gmx_resid`
        #: to transform a crystal resid to a gromacs resid.
        #:
        #: The attribute may be changed directly after init.
        self.offset = offset

        #: Auto-labelled groups use this counter.
        self.command_counter = 0

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
        self.names = self.indexfiles.keys()

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
              the individual groups: "|" (OR) or "&" (AND) ["|"]
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
        """
        if name_all is None:
            name_all = self.name_all
        if name_all is None:
            name_all = operation.join(self.indexfiles)
        if not operation in ('|', '&'):
            raise ValueError("Illegal operation %r, only '|' (OR) and '&' (AND) allowed." % 
                             operation)
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

        if len(self.selections) == 1:
            # no need to combine selections
            try:
                cmd = ['', 'q']
                rc,out,err = self.make_ndx(n=ndxfiles, o=out_ndx, input=cmd)
                if self._is_empty_group(out):
                    warnings.warn("No atoms found for %(cmd)r" % vars(), 
                                  category=BadParameterWarning)
            finally:
                if defaultgroups:
                    utilities.unlink_gmx(default_ndx)
        else:
            # multiple selections: combine them and name them
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
                name_cmd = ["name %d %s" % (last['nr'], name_all), 
                            'q']
                rc,out,err = self.make_ndx(n=tmp_ndx, o=out_ndx, input=name_cmd)
                # For debugging, look at out and err or set stdout=True, stderr=True
                # TODO: check out if at least 1 atom selected
                ##print "DEBUG: combine()"
                ##print out
            finally:
                utilities.unlink_gmx(tmp_ndx)
                if defaultgroups:
                    utilities.unlink_gmx(default_ndx)
        
        return name_all, out_ndx

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

        if selection.startswith('@'):
            process = self._process_command
            selection = selection[1:]
        else:
            process = self._process_residue
        return process(selection, name)

    def _process_command(self, command, name=None):
        """Process ``make_ndx`` command and  return name and temp index file."""
    
        self.command_counter += 1    
        if name is None:
            name = "CMD%03d" % self.command_counter
    
        # Need to build it with two make_ndx calls because I cannot reliably
        # name the new group without knowing its number.
        try:
            fd, tmp_ndx = tempfile.mkstemp(suffix='.ndx', prefix='tmp_'+name+'__')
            cmd = [command, '', 'q']   # empty command '' necessary to get list
            rc,out,err = self.make_ndx(o=tmp_ndx, input=cmd)
            self.check_output(out, "Not atoms found for selection %(command)r." % vars())
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
                        (\S\S\S)                   # 3-letter residue name
                 )
                 (?P<resid>\d+)                    # resid
                 (:                                # separator ':'
                   (?P<atom>\w+)                   # atom name
                 )?                                # possibly one 
            """, re.VERBOSE)
    
    def _process_residue(self, selection, name=None):
        """Process residue/atom selection and return name and temp index file."""

        if name is None:
            name = selection.replace(':', '_')
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

    def check_output(self, make_ndx_output, message=None):
        """Simple tests to flag problems with a ``make_ndx`` run."""
        if message is None:
            message = ""
        else:
            message = '\n' + message
        def format(output, w=60):
            hrule = "====[ GromacsError (diagnostic output) ]".ljust(w,"=")
            return hrule + '\n' + output + hrule

        rc = True
        if self._is_empty_group(make_ndx_output):
            warnings.warn("Selection produced empty group.%(message)s" 
                          % vars(), category=GromacsWarning)
            rc = False
        if self._has_syntax_error(make_ndx_output):
            rc = False
            out_formatted = format(make_ndx_output)
            raise GromacsError("make_ndx encountered a Syntax Error, "
                               "%(message)s\noutput:\n%(out_formatted)s" % vars())
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

    1. Write compact xtc and tpr with water removed.
    """

    def __init__(self, s="topol.tpr", f="traj.xtc", n=None, dirname=os.path.curdir):
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
           
        """

        self.tpr = s
        self.xtc = f
        self.ndx = n
        self.dirname = dirname

    def strip_water(self, s=None, o=None, resn="SOL", groupname="notwater", **kwargs):
        """Write compact xtc and tpr with water (by resname) removed.

        :Keywords:
           *s*
              Name of the output tpr file; by default use the original but
              insert "nowater" before suffix.
           *o*
              Name of the output trajectory; by default use the original but
              insert "nowater" before suffix.
           *resn*
              Residue name of the water molecules; all these residues are excluded.
           *groupname*
              Name of the group that is generated by subtracting all waters
              from the system.
           *kwargs* 
              are passed on to :func:`gromacs.cbook.trj_compact` (unless the
              values have to be set to certain values such as s, f, n, o
              keywords). The *input* keyword is always mangled: Only the first
              entry (the group to centre the trajectory on) is kept, and as a
              second group (the output group) *groupname* is used.

        .. warning:: The input tpr file should *not* have *any position restraints*;
                     otherwise Gromacs will throw a hissy-fit and say 

                     *Software inconsistency error: Position restraint coordinates are
                     missing*

                     (This appears to be a bug in Gromacs 4.x.)
        """
        
        def newname(name, ext, default):
            if name is None:
                p, ext = os.path.splitext(default)
                name = self.filename(p+"_nowater", ext=ext)
            return name

        newtpr = newname(s, 'tpr', self.tpr)
        newxtc = newname(o, 'xtc', self.xtc)

        nowater_ndx = "nowater.ndx"

        _input = kwargs.get('input', ['Protein'])
        kwargs['input'] = [_input[0], groupname]  # [center group, write-out selection]
        del _input

        NOTwater = "! r %(resn)s" % vars()  # make_ndx selection ("not water residues")
        with utilities.in_dir(self.dirname):
            # make no water index
            B = IndexBuilder(struct=self.tpr, selections=['@'+NOTwater], names=[groupname],
                             ndx=self.ndx, out_ndx=nowater_ndx)
            B.combine(defaultgroups=True)

            logger.info("TPR file without water %(newtpr)r" % vars())
            gromacs.tpbconv(s=self.tpr, o=newtpr, n=nowater_ndx, input=[groupname])

            logger.info("Trajectory without water %(newxtc)r" % vars())
            kwargs['s'] = self.tpr
            kwargs['f'] = self.xtc
            kwargs['n'] = nowater_ndx
            kwargs['o'] = newxtc
            trj_compact(**kwargs)

            logger.info("pdb and gro for visualization")            
            for ext in 'pdb', 'gro':
                try:
                    # see warning in doc ... so we don't use the new xtc but the old one
                    kwargs['o'] = self.filename(newtpr, ext=ext)
                    trj_compact(dump=0, stdout=False, stderr=False, **kwargs)  # silent
                except:
                    logger.exception("Failed building the water-less %(ext)s. "
                                     "Position restraints in tpr file (see docs)?" % vars())
            
            logger.info("strip_water() complete")

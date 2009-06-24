# $Id$
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)

"""
``gromacs.cbook`` -- Gromacs Cook Book
======================================

The ``cbook`` (cook book) module contains short recipes for tasks
that are often repeated. In the simplest case this is just one of the
gromacs tools with a certain set of default command line options.

By abstracting and collecting these invocations here, errors can be
reduced and the code snippets can also serve as canonical examples for
how to do simple things.

Canned Gromacs commands
-----------------------

Simple commands with new default options so that they solve a specific
problem:

.. function:: trj_compact([s="md.tpr", f="md.xtc", o="compact.xtc"[, ...]])

   Writes an output trajectory or frame with a compact representation
   of the system centered on the protein. It centers on the group
   "Protein" and outputs the whole "System" group.

.. function:: rmsd_backbone([s="md.tpr", f="md.xtc"[, ...]])

   Computes the RMSD of the "Backbone" atoms after fitting to the
   "Backbone" (including both translation and rotation).

Processing output
-----------------

There are cases when a script has to to do different things depending
on the output from a Gromacs tool. 

For instance, a common case is to check the total charge after
grompping a tpr file. The ``grompp_qtot`` function does just that.

.. autofunction:: grompp_qtot
.. autofunction:: parse_ndxlist

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

import re
import tempfile
import shutil

import gromacs
from gromacs import GromacsError
import tools

trj_compact = tools.Trjconv(ur='compact', center=True, boxcenter='tric', pbc='mol',
                            input=('protein','system'),
                            doc="Writes a compact representation of the system centered on the protein")

rmsd_backbone = tools.G_rms(what='rmsd', fit='rot+trans',
                            input=('Backbone','Backbone'),
                            doc="Computes RMSD of backbone after fitting to the backbone.")

def grompp_qtot(*args, **kwargs):
    """Run ``gromacs.grompp`` and return the total charge of the system::

      qtot = grompp_qtot(*args, **kwargs)

    where the arguments are the ones one would pass to ``gromacs.grompp``.

    :Note:

     * The stdout output of grompp is not shown. This can make debugging
       pretty hard.  Try running the normal ``gromacs.grompp`` command and
       analyze the output if the debugging messages are not sufficient.
     * Check that qtot is correct; because the function is based on pattern 
       matching of the output it canb reak when the output format changes.
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
    return qtot


def edit_mdp(mdp, new_mdp=None, **substitutions):
    """Change values in a Gromacs mdp file::

      edit_mdp('md.mdp', new_mdp='long_md.mdp', nsteps=100000, nstxtcout=1000, lincs_iter=2)

    Parameters and values are supplied as substitutions, eg ``nsteps=1000``.
    
    By default the template mdp file is **overwritten in place**.

    :Arguments:
        mdp : filename
            filename of input (and output filename of ``new_mdp=None``)
        new_mdp : filename
            filename of alternative output mdp file [None]
        substitutions : dict
            parameter=value pairs, where parameter is defined by the Gromacs mdp file

    :Returns:    
       List of parameters that have NOT been substituted.

    :Note:
       * Dashes in Gromacs mdp parameters have to be replaced by an underscore
         when supplied as python keyword arguments (a limitation of python). For example
         the MDP syntax is  ``lincs-iter = 4`` but the corresponding  keyword would be 
         ``lincs_iter = 4``.
       * If the keyword is set as a dict key, eg ``mdp_params['lincs-iter']=4`` then one
          does not have to substitute.

    :Bugs:
       * Parameters *aa_bb* and *aa-bb* are considered the same (although this should not be 
         a problem in practice because there are no mdp parameters that only differe by a underscore).
       * This code is more compact in ``Perl`` as one can use ``s///`` operators:
         ``s/^(\s*$key\s*=\s*).*/$1${val}/``

    """
    if new_mdp is None:
        new_mdp = mdp
    params = substitutions.keys()[:]   # list will be reduced for each match

    def demangled(p):
        return p.replace('_', '[-_]')  # must catch either - or _

    patterns = dict([(parameter,
                      re.compile("""\
                       (?P<assignment>\s*%s\s*=\s*)  # everything before the value
                       (?P<value>[^;]*)              # value (stop before comment=;)
                       (?P<comment>\s*;.*)?          # optional comment           
                       """ % demangled(parameter), re.VERBOSE))
                     for parameter in substitutions])

    target = tempfile.TemporaryFile()
    with open(mdp) as src:
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
    return params

def edit_txt(filename, substitutions, newname=None):
    """Primitive top editor (sed is better...)::

        edit_txt(filename, substitutions, newname=otherfilename)

        substitutions ::= [ search_replace_tuple, ... ]
        search_replace_tuple ::= ( line_match_RE, search_RE, replacement )

        line_match_RE     regular expression that selects the line (uses match)
        search_RE         regular expression that is searched in the line
        replacement       replacement string for search_RE


    ``edit_txt()`` does pretty much what a simple::

         sed /line_match_RE/s/search_RE/replacement/

    with repeated substitution commands does.

    :Note:
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
        for line in src:
            new_line = line[:]
            for subst in _substitutions:
                m = subst['lRE'].match(line)    
                if m:              # apply substition to this line?
                    #print 'match:    '+line
                    new_line = subst['sRE'].sub(subst['repl'], line)
                    #print 'replaced: '+new_line
                    break   # only apply the first matching substitution!
            target.write(new_line)
    target.seek(0)

    with open(newname, 'w') as final:
        shutil.copyfileobj(target, final)
    target.close()

    

def parse_ndxlist(output):
    """Parse output from make_ndx to build list of index groups::

      groups = parse_ndxlist(output)

    output should be the standard output from ``make_ndx``, e.g.::

       rc,output,junk = gromacs.make_ndx(..., input=('', 'q'), stdout=False, stderr=True)

    :Returns:
       The function returns a list of dicts (``groups``) with fields

       name
           name of the groups
       nr
           number of the group (starts at 0)
       natoms
           number of atoms in the group
    """
    
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
    m = NDXLIST.search(output)
    grouplist = m.group('LIST')
    NDXGROUP = re.compile(r"""
                         \s*(?P<GROUPNUMBER>\d+)      # group number
                         \s+(?P<GROUPNAME>[^\s]+)\s*: # group name, separator ':'
                         \s*(?P<NATOMS>\d+)\satoms    # number of atoms in group
                         """, re.VERBOSE)
    groups = []
    for line in grouplist.split('\n'):
        m = NDXGROUP.match(line)
        if m:
            d = m.groupdict()
            groups.append({'name': d['GROUPNAME'],
                           'nr': int(d['GROUPNUMBER']),
                           'natoms': int(d['NATOMS'])})
    return groups


# $Id$
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)

"""\
Gromacs Cook Book (cbook)
=========================

The ``cbook `` module contains short recipes for tasks that are often
repeated. In the simplest case this is just one of the gromacs tools
with a certain set of default command line options.

By abstracting and collecting these invocations here, errors can be
reduced and the code snippets can also serve as canonical examples for
how to do simple things.
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

import re
import tempfile
import shutil

import gromacs
from gromacs import GromacsError
import tools

trj_compact = tools.Trjconv(ur='compact', center=True, boxcenter='tric', pbc='mol',
                            input=('protein','system'),
                            doc="Returns a compact representation of the system centered on the protein")

rmsd_backbone = tools.G_rms(what='rmsd', fit='rot+trans',
                            input=('Backbone','Backbone'),
                            doc="Compute RMSD of backbone after fitting to the backbone.")

def grompp_qtot(*args, **kwargs):
    """Run ``gromacs.grompp`` and return the total charge of the system.

    :Bugs:
    * The stdout output of grompp is not shown. This can make debugging
      pretty hard.
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
    """Change values in a Gromacs mdp file.

      edit_mdp('md.mdp', new_mdp='long_md.mdp', nsteps=100000, nstxtcout=1000, lincs_iter=2)

    Parameters and values are supplied as substitutions, eg nsteps=1000.
    
    By default the template mdp file is **overwritten in place**.

    :Arguments:

    mdp             filename of input (and output filename of new_mdp=None)
    new_mdp         filename of alternative output mdp file [None]
    substitutions   parameter=value pairs, where parameter is defined by Gromacs

    :Returns:
    List of parameters that have NOT been substituted.


    :Notes:

    * Dashes in Gromacs mdp parameters have to be replaced by an underscore
      when supplied as python keyword arguments (a limitation of python).
      For example
        MDP:      lincs-iter = 4
        keyword:  lincs_iter = 4


    :Bugs: 

    * Parameters *aa_bb* and *aa-bb* are considered the same (but should not be a problem).
    * This code is more compact in ``Perl`` as one can use ``s///`` operators:
          s/^(\s*$key\s*=\s*).*/$1${val}/
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
    """Primitive top editor (sed is better...).

        edit_txt(filename, substitutions, newname=otherfilename)

    substitutions ::= [ search_replace_tuple, ... ]
    search_replace_tuple ::= ( line_match_RE, search_RE, replacement )

    line_match_RE     regular expression that selects the line (uses match)
    search_RE         regular expression that is searched in the line
    replacement       replacement string for search_RE


    edit_txt() does pretty much what::
      sed /line_match_RE/s/search_RE/replacement/
    with repeated substitution commands does.

    :Notes:

    * No sanity checks are performed and the substitutions must supplied
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

    



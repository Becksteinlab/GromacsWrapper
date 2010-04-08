#!/usr/bin/env python
# $Id: cleanup.py 2459 2008-11-12 17:31:10Z oliver $

"""Remove files from a staged SGE python script: parses the SGE log
file to ssh into the node and remove the stage dir. Can also parse
some of my older Charmm/SGE log scripts that contain hostname and
scratch dir information.
"""

import sys,os
import sre

def init_pattern(key):
    """Return a SRE compiled pattern; the match can be accessed in the
    match object as
      m = P[key].match(string)
      m.group(key)
    """
    return sre.compile('^init\(\): %(key)s: *(?P<%(key)s>.*)$' % locals())

INIT_KEYS = ['hostname','stagedir','JOB_ID','JOB_NAME']
P = dict( [(key,init_pattern(key)) for key in INIT_KEYS])
Q = {'hostname':sre.compile('^host: *(?P<hostname>.*)$'),
     'stagedir':sre.compile('^\+\+ temp_dir=(?P<stagedir>/scratch/oliver/.*)$'),
     'WDIR':sre.compile('^WDIR: *(?P<WDIR>.*)$'),
    }
# note: the SRE's are matched patterns, ie anchored at beginning of line!!

def scan_log(logfile,P):
    STATUS = {'abort':sre.compile('(?P<abort>Abort|abort)'),
    }
    Vars = {}
    StatusVars = {}
    log = open(logfile,"r")
    print "== %(logfile)s ==" % locals()
    for line in log:
        l = line.strip()
        for key,pattern in P.items():
            m = pattern.match(l)
            if m:
                Vars[key] = m.group(key)
                break
        for key,pattern in STATUS.items():
            m = pattern.search(l)
            if m:
                StatusVars[key] = m.group(key)
    log.close()
    return Vars, StatusVars

def cleanup(logfile):
    Vars,Status = scan_log(logfile,P)

    if len(Vars) == 0:
        print "Trying older format"
        Vars,Status = scan_log(logfile,Q)

    if len(Vars) == 0 and len(Status) == 0:
        raise ValueError('No proper tags in '+logfile)
    # all data in Var (I hope)
    print "Recognized variables: %r" % Vars
    print "Status:               %r" % Status

    try:
        # fixing older scripts which had host: instead of hostname:
        if 'hostname' not in Vars:
            Vars['hostname'] = Vars['host']
         # fix dims scripts WDIR
        if 'stagedir' not in Vars:
            Vars['stagedir'] = Vars['WDIR']

        cmd = "ssh %(hostname)s rm -vr %(stagedir)s" % Vars
    except KeyError,errmsg:
        print "Variable not found (%s)" % errmsg
        if 'abort' in Status:
            print "Job was aborted, no cleaning up necessary except log file"
            print ">>> rm "+logfile
            os.unlink(logfile)
            return
        print "Probably nfs problem with job and it never ran --- will leave log file for inspection"
        return
        #raise
    print ">>> "+cmd
    os.system(cmd)

    print ">>> rm "+logfile
    os.unlink(logfile)

if __name__ == '__main__':
    usage = "usage: %s log.oXXXXX ...\n" % sys.argv[0] + '\n' + __doc__

    try:
        logfile = sys.argv[1]
    except IndexError:
        raise ValueError("No input file.\n"+usage)
    
    for logfile in sys.argv[1:]:
        cleanup(logfile)

# $Id: common.py 2682 2008-12-23 04:32:27Z oliver $
"""Common functions for the staging module."""

import os.path

def joindicts(d1,d2):
    """Join two dictionaries u and v.

    d = joindicts(d1,d2)

    Raises a KeyError exception if two keys are identical and have
    different values.
    (This is a primitive, slow implementation. Note that it is different from dict.update().)
    """
    joined = d1.copy()  # optimize by choosing the larger dict for joined
    for k,v in d2.items():
        if k in d1 and d1[k] != v:
            raise KeyError("Duplicate key '"+str(k)+"' with value '"+str(d1[k])
                           +"' in d1 but value '"+str(v)+"' in d2.")
        joined[k] = v
    return joined


def pathjoin(head,*components,**kwargs):
    """Join all path components but unlike os.path.join sanitized them.

    If sanitize == True then
    - never allow absolute paths in second args
    - remove all '..' and makes path relative to refdir (must provide refdir)

    Arguments:
    sanitize        True: full functionality    [True]
                    False: just use os.path.join
    refdir          check if components are really under refdir [None]
    """
    kwargs.setdefault('sanitize',True)
    kwargs.setdefault('refdir',None)

    # quick exit if path needs not be scrubbed
    if not kwargs['sanitize']:
        return os.path.join(head,*components)
    
    normalized = []
    for p in components:
        p = os.path.normpath(p)
        if p.startswith('//'):
            p = p[2:]
        elif p.startswith('/'):
            p = p[1:]        
        normalized.append(p)
    npath = os.path.join(*normalized)   # path from anything before head
    
    if kwargs['refdir']:
        refdir = kwargs['refdir']
        p_abs = os.path.abspath(npath)  # probably only works when cwd == refdir (for npath not /...)
        ref_abs = os.path.abspath(refdir)
        if p_abs.find(ref_abs) >= 0:
            # good: p is under refdir
            npath = p_abs.replace(ref_abs,'')
        else:
            # try shortening the path (but that's not crucial) but use path without ..!
            common = os.path.commonprefix([p_abs,ref_abs])  # common is shorter than ref_abs
            npath = p_abs.replace(common,'',1)              # only replace at beginning else disaster with '/'!
        if npath.startswith('/'):
            npath = npath[1:]      # just in case we could not take off anything from the top

    return os.path.normpath(os.path.join(head,npath))

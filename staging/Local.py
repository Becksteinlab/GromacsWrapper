# $Id: Local.py 2050 2008-07-23 17:37:12Z oliver $
"""
:mod:`staging.Local` --- staging class for running local jobs
=============================================================

Ersatz framework for running a staged script without actually doing any staging.

Simply replace ::

   from staging.SunGridEngine import Job

with ::

   from staging.Local import Job

in the python run script (see :mod:`staging.SunGridEngine` for an
example script).

Description of the :class:`Job` class
-------------------------------------

.. autoclass:: Job
   :members:

"""
from staging.common import joindicts

class NullClass(object):
    # H. Krekel 2002
    # http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/119222
    def __getattr__(self,name): return lambda *x,**y: None

class Job(NullClass):
    """Job class that doesn't do anything but provides parameters as the 'real' classes do.

    job = Job(inputfiles=<dict>,outputfiles=<dict>,variables=<dict>,startdir=<PWD>)
    
    """

    def __init__(self,*args,**kwargs):
        super(Job,self).__init__(*args,**kwargs)
        self.input = kwargs.setdefault('inputfiles',{})
        self.output = kwargs.setdefault('outputfiles',{})
        self.filenames = joindicts(self.input,self.output)
        self.variables = kwargs.setdefault('variables',{})
        kwargs.setdefault('startdir',None)
        kwargs.setdefault('stagedir',None)
        self.queuingsystem = 'Local'

    def __repr__(self):
        return '<staging.Local.Job>'

    def save(self,filename):
        """Save the Job() as a pickled file.

        Restore with

           import staging.SunGridengine
           import cPickle
           job = cPickle.load(open(<filename>,'r'))
           
        """
        import cPickle
        cPickle.dump(self,open(filename,'wb'),cPickle.HIGHEST_PROTOCOL)
        

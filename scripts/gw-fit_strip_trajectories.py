#!/usr/bin/env python
usage = """usage: %prog [options] RUN_ID [RUN_ID ...]

Process the trajectories associated with RUN_ID, which is something
like MD_001, MD_002, etc --- essentially the directory name. All files
are assumed to be named PREFIX.* and live in directories *basedir*/RUN_ID.

Trajectories are generated for timesteps dt 1 ps, 100 ps, and 1000 ps.
"""


import gromacs
import gromacs.analysis
from gromacs.analysis import Simulation
from gromacs.analysis.plugins import StripWater
import os.path

import logging
logger = logging.getLogger("gromacs.app")

def MySimulation(identifier, **kwargs):
    basedir = kwargs.pop("basedir", os.path.curdir)
    prefix = kwargs.pop("prefix", "md")
    def F(ext, identifier=identifier, prefix=prefix):
        return os.path.join(basedir, str(identifier), 
                            "{0}.{1}".format(prefix, ext))
    dt = kwargs.pop("dt", [1, 100, 1000])
    kwargs['plugins'] = [StripWater(dt=dt, compact=True, fit="all"), ]
    return Simulation(tpr=F('tpr'), xtc=F('xtc'), edr=F('edr'), **kwargs)


if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser(usage=usage)
    parser.add_option("-B", "--basedir", dest="basedir", default=os.curdir,
                      metavar="DIR",
                      help="set basedir to DIR [%default]")
    parser.add_option("-p", "--prefix", dest="prefix", default="md",
                      metavar="PREFIX",
                      help="filenames are constructed from the default "
                      "prefix as PREFIX.xtc, PREFIX.tpr, etc (equivalent "
                      "to the -deffnm in Gromacs [%default]")
    parser.add_option("-n", "--ndx", dest="ndx", metavar="FILE",
                      default=None, 
                      help="custom index file that contains GROUP to "
                      "center on [%default]")
    parser.add_option("-g", "--center-group", dest="group", metavar="GROUP",
                      default="protein",
                      help="center trajectory on custom index group GROUP. "
                      "For special index groups also provide the index file "
                      "with the --ndx option [%default]")
    parser.add_option("-G", "--fit-group", dest="fitgroup", metavar="GROUP",
		      default="protein",
                      help="fit trajectory to index group but note that "
                      "this group MUST be a Gromacs-generated group such as 'Protein', 'DNA', 'backbone'; it cannot be one "
                      "from a custom index file (the old index file becomes invalid after "
                      "stripping of water and there is no way to regenerate a custom group) [%default]")
    parser.add_option("--force", dest="force", action="store_true",
                      default=False,
                      help="always regenerate trajectories, even if they already "
                      "exist (note that Gromacs will leave behind backup trajectories)")
                      
    opts, args = parser.parse_args()

    gromacs.start_logging()

    logger.info("Constructing filenames from prefix %(prefix)r", vars(opts))
    if len(args) == 0:
        logger.warn("No run identifiers were supplid, doing nuffing")

    for identifier in args:
        logger.info("Processing %(identifier)r...", vars())
        S = MySimulation(identifier, ndx=opts.ndx,
                         prefix=opts.prefix, basedir=opts.basedir)
        S.run('StripWater', centergroup=opts.group, fitgroup=opts.fitgroup, force=opts.force)
        logger.info("Completed %(identifier)r", vars())

    gromacs.stop_logging()

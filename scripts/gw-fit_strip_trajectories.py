#!/usr/bin/env python
usage = """usage: %prog [options] RUN_ID [RUN_ID ...]

Process the trajectories associated with RUN_ID, which is something
like MD_001, MD_002, etc --- essentially the directory name. All files
are assumed to be named md.* and live in directories *basedir*/RUN_ID.

"""
# not used
JOBNAMES = 'MD_001', 'MD_002', 'MD_003', 'MD_004', 'MD_005'

import gromacs
import gromacs.analysis
from gromacs.analysis import Simulation
from gromacs.analysis.plugins import StripWater
import os.path

import logging
logger = logging.getLogger("gromacs.app")


def MySimulation(identifier, **kwargs):
    basedir = kwargs.pop("basedir", os.path.curdir)
    def F(ext, identifier=identifier):
        return os.path.join(basedir, str(identifier), "md.%s"%ext)
    dt = kwargs.pop("dt", [1, 100, 1000])
    kwargs['plugins'] = [StripWater(dt=dt, compact=True, fit="all"), ]
    return Simulation(tpr=F('tpr'), xtc=F('xtc'), edr=F('edr'), **kwargs)


if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser(usage=usage)
    parser.add_option("-B", "--basedir", dest="basedir", default=os.curdir,
                      metavar="DIR",
                      help="set basedir to DIR [%default]")
    opts, args = parser.parse_args()

    gromacs.start_logging()

    if len(args) == 0:
        logger.warn("No run identifiers were supplid, doing nuffing")

    for identifier in args:
        logger.info("Processing %(identifier)r...", vars())
        S = MySimulation(identifier, basedir=opts.basedir)
        S.run('StripWater')
        logger.info("Completed %(identifier)r", vars())

    gromacs.stop_logging()

#!/usr/bin/env python
usage = """usage: %prog [options] [PREFIX]

Join all files from a simulation that was run in parts and which
did use the ``-noappend`` flag. The files are presumed to be named
PREFIX.part0001.EXT PREFIX.part0002.EXT ... where 
EXT = {xtc, trr, edr, log}.

The script also tries to properly name the final output frame but it
relies on the following conventions: (1) PREFIX.pdb or PREFIX.gro is
always taken as the last frame of a simulation. (2) Otherwise the
pdb/grofile with the highest parts number is taken as the final frame.

Processed files are stored in directory DIR/parts and must be
removed manually by the user; this is to avoid automatic deletion
of valuable trajectories in case something goes wrong.

The full trajectories are written to directory DIR/full.

The default for PREFIX is "md". 

.. SeeAlso:: uses :func:`gromacs.cbook.cat`, which allows access
             to more options such as specifying all directories.
"""


import gromacs.cbook
import os.path

import logging
logger = logging.getLogger("gromacs.app")



if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser(usage=usage)
    parser.add_option("-B", "--basedir", dest="basedir", default=os.curdir,
                      metavar="DIR",
                      help="find trajectories in DIR [%default]")
    opts, args = parser.parse_args()

    fulldir = "full"
    partsdir = "parts"
    fulldir_path = os.path.join(opts.basedir, fulldir)
    partsdir_path = os.path.join(opts.basedir, partsdir)

    gromacs.start_logging()

    if len(args) == 0:
        prefix = "md"
    else:
        prefix = args[0]

    logger.info("Processing PREFIX=%(prefix)r...", vars())
    try:
        gromacs.cbook.cat(prefix=prefix, dirname=opts.basedir, 
                          partsdir=partsdir, fulldir=fulldir,
                          resolve_multi="guess")
    except:
        logger.fatal("Something went wrong during joining (see below)")
        logger.fatal("To recover, manually move the processed parts from %r back to %r", partsdir_path, opts.basedir)
        logger.fatal("It is also recommended to delete %(fulldir_path)r and start from the beginning", vars())
        logger.exception("See stacktrace for details")
        raise

    logger.info("Joined parts for %(prefix)r in %(fulldir_path)r", vars())
    logger.info("Manually remove parts in %(partsdir_path)r", vars())

    gromacs.stop_logging()

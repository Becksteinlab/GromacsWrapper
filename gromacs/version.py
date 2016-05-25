# GromacsWrapper
# Copyright (c) 2009-2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.


#: Package version; this is the only place where it is set.
VERSION = 0,5,0
#: Set to ``True`` for a release. If set to ``False`` then the patch level
#: will have the suffix "-dev".
RELEASE = True
if not RELEASE:
    VERSION = VERSION[:2] + (str(VERSION[2]) + '-dev',)

def get_version():
    """Return current package version as a string."""
    return ".".join(map(str,VERSION))

def get_version_tuple():
    """Return current package version as a tuple (*MAJOR*, *MINOR*, *PATCHLEVEL*)."""
    return tuple(map(str,VERSION))

#!/usr/bin/env python
# $URL$
# $Rev$
# David Jones.
# Copyright 2008 Ravenbrook Limited.

"""Script to fetch the inputs required for the GISTEMP program.  The
inputs on are documented in the gistemp.txt file:
http://clearclimatecode.org/master/test/GISTEMP/gistemp.txt

Invoke without arguments.
"""

# http://www.python.org/doc/2.4.4/lib/module-getopt.html
import getopt
# http://www.python.org/doc/2.4.4/lib/module-sys.html
import sys
# http://www.python.org/doc/2.4.4/lib/module-urllib.html
import urllib

def doit(prefix='input/', output=sys.stdout):
    """Single function to download everything."""

    # See
    # http://groups.google.com/group/ccc-gistemp-discuss/web/compiling-gistemp-source?hl=en
    # (But station_inventory.Z corrected to station.inventory.Z)

    noaa = """
    ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/v2/v2.mean.Z
    ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/v2/v2.temperature.inv
    ftp://ftp.ncdc.noaa.gov/pub/data/ushcn/hcn_doe_mean_data.Z
    ftp://ftp.ncdc.noaa.gov/pub/data/ushcn/station.inventory.Z
    """.split()

    all = noaa

    output = sys.stdout

    # Attempt to create the directories required for prefix
    from os import makedirs
    try:
        makedirs(prefix)
    except OSError:
        # Expected if the directories already exist
        pass

    for url in all:
        def hook(n, bs, ts):
            output.write("\r%s %d" % (url, n*bs))
            output.flush()
        # Make a local filename
        local = prefix + url.split('/')[-1]
        urllib.urlretrieve(url, local, hook)
        output.write('\n')

# Guido's main, http://www.artima.com/weblogs/viewpost.jsp?thread=4829
class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "", ["help"])
            for o,a in opts:
                if o in ('--help',):
                    print __doc__
                    return 0
        except getopt.error, msg:
             raise Usage(msg)
        doit()
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, "for help use --help"
        return 2
    return 0

if __name__ == "__main__":
    sys.exit(main())

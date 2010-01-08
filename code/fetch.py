#!/usr/bin/env python
# $URL$
# $Rev$
#
# fetch.py
#
# David Jones, Ravenbrook Limited, 2010-01-07
# Copyright (C) 2008-2010 Ravenbrook Limited.

"""
fetch.py [--help] [input-files] ...

Script to fetch (download from the internet) the inputs required for
the GISTEMP program.

Any input files passed as arguments will be fetched from their locations
on the internet.  If no arguments are supplied, nothing will be fetched.
"""

# http://www.python.org/doc/2.4.4/lib/module-getopt.html
import getopt
# http://www.python.org/doc/2.4.4/lib/module-sys.html
import sys
# http://www.python.org/doc/2.4.4/lib/module-urllib.html
import urllib

def fetch(files, prefix='input/', output=sys.stdout):
    """Download a bunch of files.  *files* is a list of the files to
    download (just the basename part, the bit after the last directory
    separator, with or without the compression suffix).  *prefix* (not
    normally changed from the default) is the location on the local
    file system where the files will be downloaded.  *output* (not
    normally changed from the default) is where progress messages
    appear.
    """

    # See
    # http://groups.google.com/group/ccc-gistemp-discuss/web/compiling-gistemp-source?hl=en
    # (But station_inventory.Z corrected to station.inventory.Z)

    # There appears to be an HTTP server for the NOAA data:
    # http://www1.ncdc.noaa.gov/pub/data/ghcn/v2/
    # drj can't find any documentation that says this is the right place to
    # use.  HTTP would be better because we can find out the length of
    # the file so we can show progress better (and possibly diagnose /
    # restart interrupted downloads).

    # We are moving from USHCNv1 (v2.mean.Z) to USHCNv2
    # (9641C_200907_F52.avg.gz).  It does no harm to remember how to
    # fetch the older file.
    noaa = """
    ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/v2/v2.mean.Z
    ftp://ftp.ncdc.noaa.gov/pub/data/ushcn/v2/monthly/9641C_200907_F52.avg.gz
    ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/v2/v2.temperature.inv
    ftp://ftp.ncdc.noaa.gov/pub/data/ushcn/hcn_doe_mean_data.Z
    ftp://ftp.ncdc.noaa.gov/pub/data/ushcn/station.inventory.Z
    """.split()

    giss = """
    ftp://data.giss.nasa.gov/pub/gistemp/SBBX.HadR2
    """.split()

    all = noaa + giss

    # *name* is a dictionary that maps from short name, with and
    # without compression suffix, to the URL.
    name = {}
    for url in all:
        file = url.split('/')[-1]
        split = file.split('.')
        if split[-1].lower() in ['z','gz']:
            name['.'.join(split[:-1])] = (file, url)
        name[file] = (file, url)

    for n in files:
        if n not in name:
            raise Error("Don't know how to fetch %s" % n)

    # Attempt to create the directories required for *prefix*.
    from os import makedirs
    try:
        makedirs(prefix)
    except OSError:
        # Expected if the directories already exist.
        pass

    for file in files:
        assert file in name
        def hook(n, bs, ts):
            """Hook function for urllib.urlretrieve.  Implements a
            progress indicator.
            """

            got = n*bs
            if ts < 0:
                outof = ''
            else:
                # On the last block n*bs can exceed ts, so we clamp it
                # to avoid awkward questions.
                got = min(got, ts)
                outof = '/%d [%d%%]' % (ts, 100*got//ts)
            output.write("\r  %d%s" % (got, outof))
            output.flush()

        # Make a local filename
        (filename, url) = name[file]
        local = prefix + filename
        output.write(url + '\n')
        urllib.urlretrieve(url, local, hook)
        output.write('\n')
        output.flush()

class Error(Exception):
    """Some sort of problem with fetch."""

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
        fetch(args)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, "for help use --help"
        return 2
    return 0

if __name__ == "__main__":
    sys.exit(main())

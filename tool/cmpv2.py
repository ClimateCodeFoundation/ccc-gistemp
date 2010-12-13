#!/usr/bin/env python
# $URL$
# $Rev$
#
# cmpv2.py
#
# David Jones, Climate Code Foundation, 2010-12-13

"""
Tool to compare two files in GHCN v2 format.
"""

import gio

def compare(path1, path2, out):
    # *reader*, *station* and so are each pairs (one element for each
    # input file).
    reader = [gio.GHCNV2Reader(path=p) for p in [path1, path2]]
    station = map(list, reader)
    print >> out, "Number of stations: %d :: %d" % tuple(
      len(l) for l in station)
    uid = [set(s.uid for s in l) for l in station]
    d = uid[0] - uid[1]
    print >> out, "Stations only in %s (%d)" % (path1, len(d))
    if d:
        d = list(d)[:10]
        print >> out, ' '.join(d)
    del d
    e = uid[1] - uid[0]
    print >> out, "Stations only in %s (%d)" % (path2, len(e))
    if e:
        e = list(e)[:10]
        print >> out, ' '.join(e)
    del e


def main(argv=None):
    import sys
    if argv is None:
        argv = sys.argv

    arg = argv[1:]
    compare(*arg[0:2], out=sys.stdout)

if __name__ == '__main__':
    main()

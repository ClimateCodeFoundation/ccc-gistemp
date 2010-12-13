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
    note10(out, d)
    del d
    e = uid[1] - uid[0]
    print >> out, "Stations only in %s (%d)" % (path2, len(e))
    note10(out, e)
    del e
    common = uid[0] & uid[1]
    stationdict = [dict((st.uid,st) for st in l) for l in station]
    note_diff_years(out, common, stationdict, [path1, path2])

def note_diff_years(out, uids, station, path):
    """For the collection of uids in *uids*, note which stations in the
    pair of dicts *station* have different years."""

    count = 0
    for uid in uids:
        pair = [d[uid] for d in station]
        year = [station_years(st) for st in pair]
        if year[0] != year[1]:
            count += 1
            if count > 10:
                print >> out, "... and more"
                return
            d = year[0] - year[1]
            if d:
                print >> out, "Years only in %s %s (%d)" % (
                  uid, path[0], len(d))
                note10(out, map(str, d))
            e = year[1] - year[0]
            if e:
                print >> out, "Years only in %s %s (%d)" % (
                  uid, path[1], len(e))
                note10(out, map(str, e))

def station_years(station):
    """Return a collection of the years for which the station has
    data."""

    return set(y//12 for y in station.asdict())
        

def note10(out, s):
    """if *s* has any elements, print out up to 10 as a string on
    *out*."""
    if not s:
        return
    l = list(s)[:10]
    suffix = ''
    if len(l) < len(s):
        suffix = ' ...'
    print >> out, ' '.join(l) + suffix


def main(argv=None):
    import sys
    if argv is None:
        argv = sys.argv

    arg = argv[1:]
    compare(*arg[0:2], out=sys.stdout)

if __name__ == '__main__':
    main()

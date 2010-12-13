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
    note_diff_months(out, common, stationdict, [path1, path2])

def note_diff_years(out, uids, station, path):
    """For the collection of uids in *uids*, note which stations in the
    pair of dicts *station* have different years."""

    count = 0
    for uid in uids:
        pair = [d[uid] for d in station]
        year = [station_years(st) for st in pair]
        count = note_diffs(out, uid, year[0], year[1], path, "Years", count)
        if count > 10:
            return

def note_diffs(out, uid, a, b, path, noun, count):
    """Note differences in two sets of *noun*."""
    if a != b:
        count += 1
        if count > 10:
            print >> out, "... and more"
            return count
        d = a - b
        if d:
            print >> out, "%s only in %s %s (%d)" % (
              noun, uid, path[0], len(d))
            note10(out, map(str, d))
        e = b - a
        if e:
            print >> out, "%s only in %s %s (%d)" % (
              noun, uid, path[1], len(e))
            note10(out, map(str, e))
    return count

def station_years(station):
    """Return a collection of the years for which the station has
    data."""

    return set(y//100 for y in station.asdict())

def note_diff_months(out, uids, station, path):
    """Arguments same as *note_diff_years*; notes which stations have
    different months (assuming that periods of different years have
    already been noted).
    """

    count = 0
    for uid in uids:
        pair = [d[uid] for d in station]
        year = [station_years(st) for st in pair]
        common_years = set(year[0]) & set(year[1])
        month = [set(k for k in st.asdict() if k//100 in common_years)
          for st in pair]
        count = note_diffs(out, uid, month[0], month[1], path, "Months", count)
        if count > 10:
            return
        

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

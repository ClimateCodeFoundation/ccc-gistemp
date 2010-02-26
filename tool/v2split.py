#!/usr/bin/env python
# $URL$
# $Rev$
#
# v2split.py
#
# David Jones, Ravenbrook Limited, 2010-02-26

"""
python ghcn_split.py YYYY

Splits a GHCN v2.mean file, on stdin, into two files: v2.mean-preYYYY
v2.mean-postYYYY.  The split is made on the basis of which stations are
reporting in the year YYYY or a more recent year.  v2.mean-postYYYY will
contain records for all the stations that have a record in the year YYYY
or a more recent year; v2.mean-preYYYY will contain the records for all
the other stations.
"""

def split(inp, out, splitat):
    """Input flle: *inp*;
    Output files: *out* (a pair);
    The year used to split the stations: *splitat*.
    """

    import itertools

    def id11(line):
        """The 11-digit station identifier for a v2.mean record."""
        return line[:11]

    for stationid,lines in itertools.groupby(inp, id11):
        lines = list(lines)
        # Gather the set of years for which there are records (across
        # all duplicates for a single station).
        years = set(int(line[12:16]) for line in lines)
        if max(years) >= splitat:
            out[1].writelines(lines)
        else:
            out[0].writelines(lines)

def main(argv=None):
    import sys
    if argv is None:
        argv = sys.argv

    year = int(argv[1])
    out = [open('v2.mean-pre%d' % year, 'w'),
           open('v2.mean-post%d' % year, 'w')]

    return split(sys.stdin, out, year)

if __name__ == '__main__':
    main()

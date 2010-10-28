#!/usr/bin/env python
# $URL$
# $Rev$
#
# popchart.py
#
# David Jones, Climate Code Foundation, 2010-10-28
#
# REFERENCES
#
# P. D. Jones and A. Moberg;
# "Hemispheric and Large-Scale Surface Air Temperature Variations: An
# Extensive Revision and an Update to 2001"; Journal of Climate; 2003.

"""Tool to draw a googlechart of year versus number of data in that
year.  See Jones and Moberg 2003 Figure 1 for an example."""

import itertools
import urllib

prefix = 'http://chart.apis.google.com/chart'

def popchart(inp, out):
    """Output googlechart URL on *out*."""

    counts = count(inp)
    miny = min(counts)
    maxy = max(counts)
    most = max(counts.values())
    l = [counts.get(y, 0) for y in range(miny,maxy+1)]
    data = ','.join(map(str, l))
    xlimit = limit(most)

    d = dict(cht='bvg',
      chd='t:'+data,
      chxt='x,y,x',
      chxl='2:|year+(CE)',
      chxp='2,50',
      chs='440x330',
      chds="1,%d" % xlimit,
      chbh='a,0,2',
      chco='6611cc',
      chxr='0,%d,%d,10|1,0,%d' % (miny,maxy,xlimit),
    )

    url = prefix + '?' + '&'.join(map('='.join, d.items()))
    out.write(url+'\n')

def limit(x):
    """Calculate a limit (for the scale of the y-axis).  The result
    is either 1, 2, or 5 times some power of ten.  And is the smallest
    such number that is >= x.
    """

    if x < 2:
        return 1
    x -= 1
    s = str(x)
    if s[0] < '2':
        return 2*10**len(s[1:])
    if s[0] < '5':
        return 5*10**len(s[1:])
    return 10**len(s)

def count(inp):
    """*inp* is a file in v2.mean format.  Counts the number of rows in
    each year."""

    def getyear(row):
        return row[12:16]
    res = dict((int(year),len(list(rows))) for year,rows in
      itertools.groupby(sorted(inp, key=getyear), getyear))
    return res

def main(argv=None):
    import sys
    if argv is None:
        argv = sys.argv

    arg = argv[1:]
    if arg:
        f = open(arg[0])
    else:
        f = sys.stdin
    popchart(f, sys.stdout)

if __name__ == '__main__':
    main()

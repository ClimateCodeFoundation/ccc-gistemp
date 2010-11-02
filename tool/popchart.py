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

def popchart(inps, out):
    """Output googlechart URL on *out*. *inps* is a list of files in
    v2.mean format."""

    # There is one "count" dict for each input, the dict maps from year
    # to count of rows for that year.
    counts = [count(inp) for inp in inps]
    minyear = min(min(c) for c in counts)
    maxyear = max(max(c) for c in counts)
    most = max(max(c.values()) for c in counts)
    yscale = reasonable_scale(most)

    # Convert the counts to simple sequences (each sequence starting with
    # *minyear*.
    seqs = [[c.get(y, 0) for y in range(minyear,maxyear+1)]
      for c in counts]

    data = '|'.join(','.join(map(str, l)) for l in seqs)

    chdl=None
    try:
        chdl = '|'.join(inp.name for inp in inps)
    except:
        pass

    d = dict(cht='bvg',
      chd='t:'+data,
      chdl=chdl,
      chxt='x,y,x',
      chxl='2:|year+(CE)',
      chxp='2,50',
      chs='440x330',
      chds="1,%d" % yscale,
      chbh='a,0,2',
      chco='6611cc',
      chxr='0,%d,%d,10|1,0,%d' % (minyear,maxyear,yscale),
    )

    if not d['chdl']:
        del d['chdl']

    url = prefix + '?' + '&'.join(map('='.join, d.items()))
    out.write(url+'\n')

def reasonable_scale(x):
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
    each year.  The result is a dict that maps from year (a number) to
    count (also a number)."""

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
        fs = map(open, arg)
    else:
        fs = [sys.stdin]
    popchart(fs, sys.stdout)

if __name__ == '__main__':
    main()

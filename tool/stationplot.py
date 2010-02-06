#!/usr/bin/env python
# $URL$
# $Rev$
#
# stationplot.py
#
# David Jones, Clear Climate Code, 2010-03-04

"""
Usage: python stationplot.py [-d v2.mean] [-o file.svg] station-id

Tool to plot the records for a station.  Stations have an 11-digit
identifier, and in the GHCN file v2.mean they have a single digit added, the
duplicate marker, to form a 12-digit record identifier.  One station may
have several "duplicate" records associated with it.

Specifying an 11-digit identifier will plot all the duplicates
associated with that station; a 12-digit indentifier will plot only a
single record.

Normally the output is an SVG document written to the file plot.svg.
The -o option can be used to change this ("-o -" specifies stdout).

Normally the input is the GHCN dataset input/v2.mean; the -d option can
be used to change this ("-d -" specifies stdin).
"""

import os
import sys

# :todo: yukh!
# Clear Climate Code
sys.path.append(os.path.join(os.getcwd(),'code'))

# The hex colours come from
# http://www.personal.psu.edu/cab38/ColorBrewer/ColorBrewer.html
# (and reordered)
colour_list = """
blue
deeppink
green
orange
olive
navy
aqua
fuchsia
gray
#1f77b4
#33a02c
#fb9a99
#cab2d6
#e31a1c
#fdbf6f
#ff7f00
#6a3d9a
#b2df8a
#a6cee3
""".split()

def aplot(rows):
    """Take the lines in `rows` and splits the data into contiguous
    lists suitable for plotting.  A sequence of lists is returned, each
    list being a list of (year, datum) pairs.  Curiously year will be a
    fractional year."""

    from itertools import groupby
    import struct

    BAD = -9999

    def point(year, month, datum):
        """Convert year, month, datum to a coordinate."""
        # Put each data point on the centre of its month, assuming
        # all months are of equal width.
        return (year + (month+0.5)/12.0, int(datum))

    def asstream():
        # previous year
        prev = None
        for row in rows:
            year = float(row[12:16])
            if prev:
                assert year > prev
                if prev+1 < year:
                    # Entire missing year in record, we only need one
                    # BAD datam to force a break.
                    yield point(prev, 0, BAD)
            prev = year
            for i,datum in enumerate(struct.unpack('5s'*12, row[16:-1])):
                yield point(year, i, datum)
    for isbad,block in groupby(asstream(), lambda x: x[1] == BAD):
        if isbad:
            continue
        yield list(block)

def plot(arg, inp, out):
    """Read data from `inp` and create a plot of the stations specified
    in the list `arg`.  Plot is written to `out`.
    """

    import struct

    BAD = -9999

    table = asdict(arg, inp)
    minyear = 9999
    maxyear = -9999
    highest = -9999
    lowest = 9999
    for id12,lines in table.items():
        for row in lines:
            year = int(row[12:16])
            minyear = min(minyear, year)
            maxyear = max(maxyear, year)
            data = struct.unpack('5s'*12, row[16:-1])
            for datum in map(int, data):
                if datum == BAD:
                    continue
                highest = max(highest, datum)
                lowest = min(lowest, datum)
    # The data should be such that a station cannot have entirely
    # invalid data.  At least one year should have at least one valid
    # datum.
    assert highest > -9999
    assert lowest < 9999
    highest /= 10.0
    lowest /= 10.0

    limyear = maxyear + 1
    # Bounds of the box that displays data.  In SVG viewBox format.
    databox = (minyear, lowest, limyear-minyear, highest-lowest)
    plotwidth = databox[2]
    plotheight = databox[3]

    out.write("""<svg width='1000px' height='750px' viewBox='0 0 %d %d'
      xmlns="http://www.w3.org/2000/svg" version="1.1">\n""" %
      (plotwidth+8, plotheight+10))

    # Style
    out.write("""<defs>
  <style type="text/css">
    path { stroke-width: 0.1; fill: none }
    path.singleton { stroke-width: 0.2; stroke-linecap: round }
    g#axes path { stroke-width:0.2; fill:none; stroke: black }
    g#axes text { fill: black; font-family: Verdana }
""")
    assert len(table) <= len(colour_list)
    for id12,colour in zip(table, colour_list):
        out.write("    g#record%s { stroke: %s }\n" % (id12, colour))
    out.write("  </style>\n</defs>\n")

    # push chart down and right to give a bit of a border
    out.write("<g transform='translate(4,4)'>\n")
    # Transform so that (0,0) on chart is lower left
    out.write("<g transform='translate(0, %.1f)'>\n" % (databox[3]))
    # In this section 0,0 should coincide with bottom left of chart, but
    # oriented as per SVG default.  +ve y is down.
    out.write("<g id='axes'>\n")
    w = limyear - minyear
    # Ticks.
    s = (-minyear)%10
    # Where we want ticks, in years offset from the earliest year.
    # We have ticks every decade.
    tickat = range(s, w+1, 10)
    out.write("  <path d='" +
      ''.join(map(lambda x: 'M%d 0l0 2' % x, tickat)) +
      "' />\n")
    # Labels.
    # Font size.  Couldn't get this to work in the style element, so it's
    # here as an attribute on each <text>
    fs = 2
    for x in tickat:
        out.write("  <text text-anchor='middle'"
          " font-size='%.1f' x='%d' y='4'>%d</text>\n" %
          (fs, x, minyear+x))
    # End of "axes" group.
    out.write("</g>\n")

    # Transform so that up (on data chart) is +ve.
    out.write("<g transform='scale(1, -1)'>\n")
    # Transform so that databox lower left ends up at 0,0
    out.write("<g transform='translate(%d, %.1f)'>\n" %
      (-databox[0], -databox[1]))
    out.write("""<rect x='%d' y='%.1f' width='%d' height='%.1f'
      stroke='pink' fill='none' opacity='0.30' />\n""" % databox)

    for id12,lines in table.items():
        out.write("<g id='record%s'>\n" % id12)
        for segment in aplot(lines):
            out.write(aspath(segment)+'\n')
        out.write("</g>\n")
    out.write("</g>\n" * 4)
    out.write("</svg>\n")

def aspath(l):
    """Encode a list of data points as an SVG path element.  The element
    is returned as a string."""

    assert len(l) > 0

    # Format an (x,y) tuple.
    def fmt(t):
        return "%.3f %.1f" % (t[0], t[1]/10.0)

    d = 'M'+fmt(l[0])+'L'+' '.join(map(fmt, l[1:]))
    decorate = ''
    if len(l) == 1:
        # For singletons we:
        # - draw a length 0 segment to force a real stroke;
        # - add a class attribute so that they can be styled with larger
        # blobs.
        assert d[-1] == 'L'
        d = d[:-1] + 'l 0 0'
        decorate = "class='singleton' "
    return "<path %sd='%s' />" % (decorate, d)
        

def asdict(arg, inp):
    """`arg` should be a list of 11-digit station identifiers or
    12-digit record identifiers.  The records from `inp` are extracted
    and returned as a dictionary (that maps identifiers to lists of
    lines).
    """

    def id12(line):
        """Return the 12-digit record identifier for a line of the
        v2.mean file."""
        return line[:12]

    # http://www.python.org/doc/2.4.4/lib/module-itertools.html
    from itertools import groupby

    table = {}
    for id12,lines in groupby(inp, id12):
        id11 = id12[:11]
        if id12 in arg or id11 in arg:
            table[id12] = list(lines)
    return table
        
def main(argv=None):
    # http://www.python.org/doc/2.4.4/lib/module-getopt.html
    import getopt
    import sys
    if argv is None:
        argv = sys.argv

    outfile = 'plot.svg'
    infile = 'input/v2.mean'
    opt,arg = getopt.getopt(argv[1:], 'o:d:')
    for k,v in opt:
        if k == '-o':
            outfile = v
        if k == '-d':
            infile = v
    if outfile == '-':
        outfile = sys.stdout
    else:
        outfile = open(outfile, 'w')
    # :todo: yukh!
    from step0 import open_or_uncompress
    if infile == '-':
        infile = sys.stdin
    else:
        infile = open_or_uncompress(infile)
    return plot(arg, inp=infile, out=outfile)

if __name__ == '__main__':
    main()

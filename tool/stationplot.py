#!/usr/bin/env python
# $URL$
# $Rev$
#
# stationplot.py
#
# David Jones, Clear Climate Code, 2010-03-04
#
# For testing purposes it might be interesting to note the following:
# 61710384000 longest timespan
# 10160400001 briefest timespan
# 22224266000 greatest temperature variation
# 50894004001 least temperature variation
# 30781001000 widest aspect ratio
# 21544218001 tallest aspect ratio

"""
Usage: python stationplot.py [-d v2.mean] [-o file.svg] [-a] station-id

Tool to plot the records for a station.  Stations have an 11-digit
identifier, and in the GHCN file v2.mean they have a single digit added, the
duplicate marker, to form a 12-digit record identifier.  One station may
have several "duplicate" records associated with it.

Anomalies can be plotted (each datum has the mean for that month
subtracted from it) by using the -a option.

Specifying an 11-digit identifier will plot all the duplicates
associated with that station; a 12-digit indentifier will plot only a
single record.

Normally the output is an SVG document written to the file plot.svg.
The -o option can be used to change where it written ("-o -" specifies stdout).

Normally the input is the GHCN dataset input/v2.mean; the -d option
specifies an alternate input file ("-d -" specifies stdin).
"""

import math
import os
import sys

# Clear Climate Code
import extend_path

# :todo: Should really import this from somewhere.  Although this BAD
# value is entirely internal to this module.
BAD = 9999

class Error(Exception):
    """Some sort of error."""

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

def aplot(series):
    """`series` is a (data,begin) pair.  Each datum is enumerated with
    its fractional year coordinate, and the entire series is split into
    contiguous chunks.  The results are intended to be suitable for
    plotting.

    A sequence of lists is returned, each list being a list of
    (year, datum) pairs.  year will be a fractional year.
    """

    from itertools import groupby

    def enum_month(data):
        """Like enumerate, but decorating each datum with a fractional
        year coordinate."""

        for m,datum in enumerate(data):
            yield (first + (m+0.5)/12.0, datum)

    data,first = series

    for isbad,block in groupby(enum_month(data), lambda x: x[1] == BAD):
        if isbad:
            continue
        yield list(block)

def plot(arg, mode, inp, out, meta):
    """Read data from `inp` and create a plot of the stations specified
    in the list `arg`.  Plot is written to `out`.  Metadata (station
    name, location) is takem from the `meta` file (usually v2.inv).
    `mode` should be 'temp' to plot temperatures, or 'anom' to plot
    monthly anomalies.
    """

    import struct

    def valid(datum):
        return datum != BAD

    table = asdict(arg, inp, mode)
    if not table:
        raise Error('No data found for %s' % (', '.join(arg)))
        
    if meta:
        meta = get_meta(table, meta)
        title = []
        for id11,d in meta.items():
            title.append('%s %+06.2f%+07.2f  %s' %
              (id11, d['lat'], d['lon'], d['name']))
    title = '\n'.join(title)

    # Calculate first and last year, and highest and lowest temperature.
    minyear = 9999
    limyear = -9999
    highest = -9999
    lowest = 9999
    for _,(data,begin) in table.items():
        minyear = min(minyear, begin)
        limyear = max(limyear, begin+len(data)//12)
        valid_data = filter(valid, data)
        ahigh = max(valid_data)
        alow = min(valid_data)
        highest = max(highest, ahigh)
        lowest = min(lowest, alow)
    # The data should be such that a station cannot have entirely
    # invalid data.  At least one year should have at least one valid
    # datum.
    assert highest > -9999
    assert lowest < 9999

    # Bounds of the box that displays data.  In SVG viewBox format.
    databox = (minyear, lowest, limyear-minyear, highest-lowest)

    # Vertical scale.  The data y-coordinate is multipled by vs to get
    # to pixels (then possibly shifted up or down).  This is thus the
    # number of pixels per degree C.
    vs = 10
    # Horizontal scale.  The data x-coordinate is multiplied by hs to
    # get to pixels (then possibly shifted left or right).  This is thus
    # the number of pixels per year.
    hs = 6

    plotwidth = databox[2] * hs

    # Bottom edge and top edge of plot area, after data has been scaled.
    # Forcing them to be integers means our plot can be aligned on
    # integer coordinates.
    ybottom = math.floor((lowest-0.05)*vs)
    ytop = math.ceil((highest+0.05)*vs)
    plotheight = ytop - ybottom

    out.write("""<svg width='%dpx' height='%dpx'
      xmlns="http://www.w3.org/2000/svg" version="1.1">\n""" %
      (plotwidth+100, plotheight+100))

    # Style
    out.write("""<defs>
  <style type="text/css">
    path { stroke-width: 0.1; fill: none }
    path.singleton { stroke-width: 0.2; stroke-linecap: round }
    g#axes path { stroke-width:1; fill:none; stroke: #888 }
    g#axes text { fill: black; font-family: Verdana }
    g#title text { fill: black; font-family: Verdana }
""")
    assert len(table) <= len(colour_list)
    for id12,colour in zip(table, colour_list):
        out.write("    g#record%s { stroke: %s }\n" % (id12, colour))
    out.write("  </style>\n</defs>\n")

    # push chart down and right to give a bit of a border
    out.write("<g transform='translate(80,80)'>\n")

    # In this section 0,0 is at top left of chart, and +ve y is down.
    if title:
        out.write("  <g id='title'>\n")
        out.write("  <text font-size='20' x='0' y='-4'>%s</text>\n" % title)
        out.write("  </g>\n")

    # Transform so that (0,0) on chart is lower left
    out.write("<g transform='translate(0, %.1f)'>\n" % plotheight)
    # In this section 0,0 should coincide with bottom left of chart, but
    # oriented as per SVG default.  +ve y is down.

    # Start of "axes" group.  In this group we are 1-1 with SVG pixels;
    # (0,0) is at plot lower left, and +ve is down.  Use vs and xs to
    # scale to/and from data coordinates.
    out.write("<g id='axes'>\n")
    w = limyear - minyear
    # Ticks on the horizontal axis.
    s = (-minyear)%10
    # Where we want ticks, in years offset from the earliest year.
    # We have ticks every decade.
    tickat = range(s, w+1, 10)
    # Amount by which tick shoots out at the bottom.
    overshoot = 16
    out.write("  <path d='" +
      ''.join(map(lambda x: 'M%d %.1fl0 %.1f' %
      (x*hs, overshoot, -(plotheight+overshoot)), tickat)) +
      "' />\n")
    # Labels.
    # Font size.  Couldn't get this to work in the style element, so it's
    # here as an attribute on each <text>
    fs = 16
    for x in tickat:
        out.write("  <text text-anchor='middle'"
          " font-size='%.1f' x='%d' y='%d'>%d</text>\n" %
          (fs, x*hs, overshoot, minyear+x))
    out.write("  <g id='vaxis' font-size='%.1f' text-anchor='end'>" % fs)
    # Ticks on the vertical axis.
    # Ticks every 5 degrees C
    every = 5*vs
    s = (-ybottom) % every
    tickat = map(lambda x:x+s, range(0, int(plotheight+1-s), every))
    out.write("  <path d='" +
      ''.join(map(lambda y: 'M0 %.1fl-8 0' % -y, tickat)) +
      "' />\n")
    for y in tickat:
        # Note: "%.0f' % 4.999 == '5'
        out.write("  <text alignment-baseline='middle'"
          " x='-8' y='%.1f'>%.0f</text>\n" %
          (-y, (y+ybottom)/float(vs)))
    # End of vertical axis group.
    out.write("</g>\n")
    # End of "axes" group.
    out.write("</g>\n")

    # Transform so that up (on data chart) is +ve.
    out.write("<g transform='scale(1, -1)'>\n")
    # Transform so that plot lower left is at (0,0):
    out.write("<g transform='translate(0,%.0f)'>\n" % -ybottom)
    # Transform by (hs,vs) so that outside this group we're in (SVG)
    # pixels.
    out.write("<g transform='scale(%.2f,%.2f)'>\n" % (hs, vs))
    # Transform so that databox left ends up at x=0
    out.write("<g transform='translate(%d,0)'>\n" % (-minyear))
    out.write("""<rect x='%d' y='%.1f' width='%d' height='%.1f'
      stroke='pink' fill='none' opacity='0.30' />\n""" % databox)

    for id12,series in table.items():
        out.write("<g id='record%s'>\n" % id12)
        for segment in aplot(series):
            out.write(aspath(segment)+'\n')
        out.write("</g>\n")
    out.write("</g>\n" * 6)
    out.write("</svg>\n")

def get_meta(l, meta):
    """For the 11-digit stations identifiers in `l`, get the metadata
    extracted from the file `meta`.  A dictionary is returned from maps from
    station id to an info dictionary.  The info dictionary has keys:
    name, lat, lon (and maybe more in future).
    """

    full = {}
    for line in meta:
        id = line[:11]
        full[id] = dict(
            name = line[12:42].strip(),
            lat = float(line[43:49]),
            lon = float(line[50:57]),
        )
    d = {}
    l = set(map(lambda x: x[:11], l))
    for id11 in l:
        if id11 in full:
            d[id11] = full[id11]
    return d

def aspath(l):
    """Encode a list of data points as an SVG path element.  The element
    is returned as a string."""

    assert len(l) > 0

    # Format an (x,y) tuple.
    def fmt(t):
        return "%.3f %.1f" % (t[0], t[1])

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

# Pasted from
# http://code.google.com/p/ccc-gistemp/source/browse/trunk/code/step1.py?r=251
# :todo: abstract properly.
def from_years(years):
    """*years* is a list of year records (lists of temperatures) that
    comprise a station's entire record.  The data are converted to a
    linear array (could be a list/tuple/array/sequence, I'm not
    saying), *series*, where series[m] gives the temperature (a
    floating point value in degrees C) for month *m*, counting from 0
    for the January of the first year with data.

    (*series*,*begin*) is returned, where *begin* is
    the first year for which there is data for the station.

    This code is also in step0.py at present, and should be shared.
    """

    begin = None
    # Previous year.
    prev = None
    series = []
    for (year, data) in years:
        if begin is None:
            begin = year
        # The sequence of years for a station record is not
        # necessarily contiguous.  For example "1486284000001988" is
        # immediately followed by "1486284000001990", missing out 1989.
        # Extend with blank data.
        while prev and prev < year-1:
            series.extend([BAD]*12)
            prev += 1
        prev = year
        series.extend(data)
    return (series, begin)
def from_lines(lines):
    """*lines* is a list of lines (strings) that comprise a station's
    entire record.  The lines are an extract from a file in the same
    format as the GHCN file v2.mean.Z.  The data are converted to a
    linear array (could be a list/tuple/array/sequence, I'm not saying),
    *series*, where series[m] gives the temperature (a floating
    point value in degrees C) for month *m*, counting from 0 for the
    January of the first year with data.

    (*series*,*begin*) is returned, where *begin* is
    the first year for which there is data for the station.

    Invalid data are marked in the input file with -9999 but are
    translated in the data arrays to BAD.
    """

    years = []
    # Year from previous line.
    prev = None
    # The previous line itself.
    prevline = None
    for line in lines:
        year = int(line[12:16])
        if prev == year:
            # There is one case where there are multiple lines for the
            # same year for a particular station.  The v2.mean input
            # file has 3 identical lines for "8009991400101971"
            if line == prevline:
                print "NOTE: repeated record found: Station %s year %s; data are identical" % (line[:12],line[12:16])
                continue
            # This is unexpected.
            assert 0, "Two lines specify different data for %s" % line[:16]
        # Check that the sequence of years increases.
        assert not prev or prev < year

        prev = year
        prevline = line
        temps = []
        for m in range(12):
            datum = int(line[16+5*m:21+5*m])
            if datum == -9999:
                datum = BAD
            else:
                # Convert to floating point and degrees C.
                datum *= 0.1
            temps.append(datum)
        years.append((year, temps))
    return from_years(years)
        

def asdict(arg, inp, mode):
    """`arg` should be a list of 11-digit station identifiers or
    12-digit record identifiers.  The records from `inp` are extracted
    and returned as a dictionary (that maps identifiers to (data,begin)
    pair).  If `mode` is 'anom' then data are converted to monthly
    anomalies.
    """

    # Clear Climate Code, tool directory
    import v2index

    v2 = v2index.File(inp)

    table = {}
    for id in arg:
        for id12,series in v2.get(id):
            data,begin = from_lines(series)
            if mode == 'anom':
                # Clear Climate Code, code directory
                from code.step1 import month_anomaly
                _,data = month_anomaly(data)
            table[id12] = (data,begin)

    return table

class Usage(Exception):
    pass
        
def main(argv=None):
    # http://www.python.org/doc/2.4.4/lib/module-getopt.html
    import getopt
    import sys
    if argv is None:
        argv = sys.argv

    try:
        outfile = 'plot.svg'
        infile = 'input/v2.mean'
        metafile = 'input/v2.inv'
        mode = 'temp'
        opt,arg = getopt.getopt(argv[1:], 'ao:d:m:')
        if not arg:
            raise Usage('At least one identifier must be supplied.')
        for k,v in opt:
            if k == '-a':
                mode = 'anom'
            if k == '-o':
                outfile = v
            if k == '-d':
                infile = v
            if k == '-m':
                metafile = v
        if outfile == '-':
            outfile = sys.stdout
        else:
            outfile = open(outfile, 'w')
        if infile == '-':
            infile = sys.stdin
        else:
            infile = open(infile)
        metafile = open(metafile)
        return plot(arg, mode=mode, inp=infile, out=outfile, meta=metafile)
    except (getopt.GetoptError, Usage), e:
        sys.stdout.write('%s\n' % str(e))
        sys.stdout.write(__doc__)

if __name__ == '__main__':
    main()

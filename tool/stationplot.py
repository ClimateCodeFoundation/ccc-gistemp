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
Usage: python stationplot.py [options] station-id

The options are:
  [-a] [-y] [-d v2.mean] [-o file.svg] [-t YYYY,YYYY] [-s 0.01] [-c config]

Tool to plot the records for a station.  Stations have an 11-digit
identifier, and in the GHCN file v2.mean they have a single digit added, the
duplicate marker, to form a 12-digit record identifier.  One station may
have several "duplicate" records associated with it.

Specifying an 11-digit identifier will plot all the duplicates
associated with that station; a 12-digit indentifier will plot only a
single record.

Anomalies can be plotted (each datum has the mean for that month
subtracted from it) by using the -a option.  The -y option will compute
monthly anomalies and then average them in blocks of 12 to form yearly
(annual) anomalies.

The -t option will restrict the time axis so that only records between
the beginning of the first year and the beginning of the second year are
displayed (in other words, the second year is excluded from the
displayed range).

Normally the output is an SVG document written to the file
station-id.svg (in other words, the first argument with ".svg" appended).
The -o option can be used to change where it written ("-o -" specifies stdout).

Normally the input is the GHCN dataset input/v2.mean; the -d option
specifies an alternate input file ("-d -" specifies stdin).

Normally the input data is expected to be temperatures in units of 0.1 C
(this is the convention used by GHCN v2).  The -s option can be used to
specify a different sized unit.

The -c option can be used to set various configuration options.  Best to
examine the source code for details.
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

class Config:
    """A record of the configuration parameters used to style the plot.
    Just a struct really."""

config = Config()
config.debug = False
config.fontsize = 16
# Amount by which tick shoots out at the bottom.
config.overshoot = 16
# Pixels per year.
config.xscale = 6
# Pixels per degree C.
config.yscale = 10
# Ticks on Y axis (in scaled data units, that is, degrees C).
config.ytick = 5
# Style of legend.  'pianola' or 'none'.
config.legend = 'pianola'

def derive_config(config):
    """Some configuration parameters are derived from others if they
    haven't been set."""

    def titlesize(c):
        c.titlesize = 1.25*c.fontsize

    d = dict(titlesize=titlesize)
    for attr,fn in d.items():
        if not hasattr(config, attr):
            fn(config)


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

def aplot(series, K):
    """`series` is a (data,begin) pair.  Each datum is enumerated with
    its fractional year coordinate, and the entire series is split into
    contiguous chunks.  The results are intended to be suitable for
    plotting.

    *K* is the number of data items per year.  This is 12 for monthly
    data; 1 for annual data.

    A sequence of lists is returned, each list being a list of
    (year, datum) pairs.  year will be a fractional year.
    """

    from itertools import groupby

    def enum_frac(data):
        """Like enumerate, but decorating each datum with a fractional
        year coordinate."""

        for m,datum in enumerate(data):
            yield (first + (m+0.5)/K, datum)

    data,first = series

    for isbad,block in groupby(enum_frac(data), lambda x: x[1] == BAD):
        if isbad:
            continue
        yield list(block)

def plot(arg, inp, out, meta, timewindow=None, mode='temp', scale=0.1):
    """Read data from `inp` and create a plot of the stations specified
    in the list `arg`.  Plot is written to `out`.  Metadata (station
    name, location) is taken from the `meta` file (usually v2.inv).
    `mode` should be 'temp' to plot temperatures, or 'anom' to plot
    monthly anomalies.  `timewindow` restricts the plot to a particular
    range of times: None means that the entire time range is plotted;
    otherwise, it should be a pair of numbers (y1,y2) and only records
    that have a time t where y1 <= t < y2 are displayed.  Normally y1
    and y2 are years in which case records from the beginning of y1 up
    to the beginning of y2 are displayed.
    """

    import itertools
    import struct

    def valid(datum):
        return datum != BAD

    datadict = asdict(arg, inp, mode, scale)
    if not datadict:
        raise Error('No data found for %s' % (', '.join(arg)))
        
    if meta:
        meta = get_meta(datadict, meta)
        title = []
        for id11,d in meta.items():
            title.append('%s %+06.2f%+07.2f  %s' %
              (id11, d['lat'], d['lon'], d['name']))
    title = '\n'.join(title)

    # Assign number of data items per year.
    global K # :todo: made not global.
    if 'annual' in mode:
        K = 1
    else:
        K = 12

    # Calculate first and last year, and highest and lowest temperature.
    minyear = 9999
    limyear = -9999
    highest = -9999
    lowest = 9999
    if timewindow:
        # :todo: make work with mode=annual
        datadict = window(datadict, timewindow)
    for _,(data,begin) in datadict.items():
        minyear = min(minyear, begin)
        limyear = max(limyear, begin + (len(data)//K))
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

    plotwidth = databox[2] * config.xscale

    # Bottom edge and top edge of plot area, after data has been scaled.
    # Forcing them to be integers means our plot can be aligned on
    # integer coordinates.
    ybottom = math.floor((lowest-0.05)*config.yscale)
    ytop = math.ceil((highest+0.05)*config.yscale)
    plotheight = ytop - ybottom
    legendh = legend_height(datadict)
    lborder = 125
    rborder = config.fontsize*2

    out.write("""<svg width='%dpx' height='%dpx'
      xmlns="http://www.w3.org/2000/svg"
      xmlns:xlink="http://www.w3.org/1999/xlink"
      version="1.1">\n""" %
      (plotwidth+lborder+rborder, plotheight+100+legendh))

    # Style
    out.write("""<defs>
  <style type="text/css">
    .debug { %s }
    path { stroke-width: 1.4; fill: none }
    path.singleton { stroke-width: 2.8; stroke-linecap: round }
    g#axes path { stroke-width:1; fill:none; stroke: #888 }
    g#legend path { stroke-width:1; fill:none }
    text { fill: black; font-family: Verdana }
""" % ('display: none', '')[config.debug])
    colours = itertools.chain(colour_list, colour_iter())
    for id12,colour in zip(datadict, colours):
        cssidescaped = cssidescape('record' + id12)
        out.write("    g.%s { stroke: %s }\n" % (cssidescaped, colour))
    out.write("  </style>\n</defs>\n")

    # Push chart down and right to give a bit of a border.
    out.write("<g transform='translate(%.1f,80)'>\n" % lborder)

    # In this section 0,0 is at top left of chart, and +ve y is down.
    if title:
        out.write("  <g id='title'>\n")
        out.write("  <text font-size='%.1f' x='0' y='-4'>%s</text>\n" %
          (config.titlesize, title))
        out.write("  </g>\n")

    # Transform so that (0,0) on chart is lower left
    out.write("<g transform='translate(0, %.1f)'>\n" % plotheight)
    # In this section 0,0 should coincide with bottom left of chart, but
    # oriented as per SVG default.  +ve y is down.  We are 1-1 with SVG
    # pixels.
    # Use yscale and xscale to scale to/and from data coordinates.

    # Colour legend.
    out.write("<g id='legend'>\n")
    render_legend(datadict, minyear, out)
    # End of legend group
    out.write("</g>\n")

    # Start of "axes" group.
    out.write("<g id='axes'>\n")
    w = limyear - minyear
    # Ticks on the horizontal axis.
    s = (-minyear)%10
    # Where we want ticks, in years offset from the earliest year.
    # We have ticks every decade.
    tickat = range(s, w+1, 10)
    out.write("  <path d='" +
      ''.join(map(lambda x: 'M%d %.1fl0 %.1f' %
      (x*config.xscale, config.overshoot, -(plotheight+config.overshoot)),
      tickat)) +
      "' />\n")
    # Horizontal labels.
    for x in tickat:
        out.write("  <text text-anchor='middle'"
          " font-size='%.1f' x='%d' y='%d'>%d</text>\n" %
          (config.fontsize, x*config.xscale, config.overshoot, minyear+x))
    # Vertical axis.
    out.write("  <g id='vaxis' font-size='%.1f' text-anchor='end'>" %
      config.fontsize)
    # Ticks on the vertical axis.
    # Ticks every *ytick* degrees C.
    every = config.ytick*config.yscale
    # Works best when *every* is an integer.
    assert int(every) == every
    every = int(every)
    s = (-ybottom) % every
    tickat = map(lambda x:x+s, range(0, int(plotheight+1-s), every))
    out.write("  <path d='" +
      ''.join(map(lambda y: 'M0 %.1fl-8 0' % -y, tickat)) +
      "' />\n")
    # *prec* gives the number of figures after the decimal point for the
    # y-axis tick label.
    prec = -int(round(math.log10(config.ytick) - 0.001))
    prec = max(0, prec)
    tickfmt = '%%.%df' % prec
    for y in tickat:
        # Note: "%.0f' % 4.999 == '5'
        out.write(("  <text alignment-baseline='middle'"
          " x='-8' y='%.1f'>" + tickfmt + "</text>\n") %
          (-y, (y+ybottom)/float(config.yscale)))
    # Vertical label.
    out.write("  <defs><path id='pvlabel' d='M-%d -20l0 -400'/></defs>\n" %
      (3.5*config.fontsize-8))
    if 'temp' in mode:
        value = 'Temperature'
    else:
        value = 'Anomaly'
    out.write("  <text text-anchor='start'>"
      "<textPath xlink:href='#pvlabel'>"
      u"%s (\N{DEGREE SIGN}C)</textPath></text>\n" % value)
    # End of vertical axis group.
    out.write("</g>\n")
    # End of "axes" group.
    out.write("</g>\n")

    # Transform so that up (on data chart) is +ve.
    out.write("<g transform='scale(1, -1)'>\n")
    # Transform so that plot lower left is at (0,0):
    out.write("<g transform='translate(0,%.0f)'>\n" % -ybottom)
    xdatabox = (0, ybottom, databox[2]*config.xscale, databox[3]*config.yscale)
    out.write("""<rect class='debug' x='%d' y='%.1f' width='%d' height='%.1f'
      stroke='pink' fill='none' opacity='0.30' />\n""" % xdatabox)

    def scale(points):
        """Take a sequence of (year,datum) pairs and scale using
        config.xscale and config.yscale to be in a coordinate system
        where (minyear,0) is at 0,0 and we're 1 to 1 with SVG pixels.
        """

        return [((x-minyear)*config.xscale, y*config.yscale) for x,y in points]

    for id12,series in datadict.items():
        out.write("<g class='record%s'>\n" % id12)
        for segment in aplot(series, K):
            out.write(aspath(scale(segment))+'\n')
        out.write("</g>\n")
    out.write("</g>\n" * 4)
    out.write("</svg>\n")

def legend_height(datadict):
    """Height of the legend.  Seems to include the x-axis label (?)"""

    if config.legend == 'pianola':
        return config.fontsize*(len(datadict)+1)
    else:
        return config.fontsize

def render_legend(datadict, minyear, out):
    """Write the SVG for the legend onto the stream *out*."""
    # Includes range indicators for each series.
    # Start of the top line of text for the legend.

    if config.legend == 'pianola':
        yleg = config.overshoot+config.fontsize
        for i,(id12,(data,begin)) in enumerate(sorted(datadict.items())):
            length = len(data)//K
            y = yleg + config.fontsize*i
            out.write("  <text alignment-baseline='middle'"
              " text-anchor='end' x='0' y='%d'>%s</text>\n" %
              (y, id12))
            out.write("  <g class='record%s'><path d='M%.1f %dl%.1f 0' /></g>\n" %
              (id12, (begin-minyear)*config.xscale, y, length*config.xscale))

def cssidescape(identifier):
    """Escape an identifier so that it is suitable for use as a CSS
    identifier.  See http://www.w3.org/TR/CSS2/syndata.html ."""

    import re

    def f(m):
       return '\\'+m.group()

    x0 = identifier[0]

    # Escape all but initial character.
    x = re.sub(r'([^a-zA-Z0-9_-])', f, identifier[1:])
    # Initial character needs escaping too.
    if x0 in '0123456789':
        x0 = '\\%02x' % ord(x0)
    elif not re.match(r'^[a-zA-Z_]', x0):
        x0 = '\\' + x0
    return x0 + x

def window(datadict, timewindow):
    """Restrict the data series in *datadict* to be between the two
    times specified in the timewindow pair.  A fresh dict is returned.
    """

    t1,t2 = timewindow
    # The window must be on a year boundary to preserve the fact that
    # data is a multiple of 12 long.
    assert int(t1) == t1
    assert int(t2) == t2
    d = {}
    for id12,(data,begin) in datadict.items():
        if t2 <= begin:
            continue
        end = begin+len(data)//12
        if end <= t1:
            continue
        if t2 < end:
            data = data[:12*(t2-end)]
        if begin < t1:
            data = data[12*(t1-begin):]
            begin = t1
        d[id12] = (data,begin)
    return d

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

def from_lines(lines, scale=0.1):
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
                datum *= scale
            temps.append(datum)
        years.append((year, temps))
    return from_years(years)
        

def asdict(arg, inp, mode, scale=0.1):
    """`arg` should be a list of 11-digit station identifiers or
    12-digit record identifiers.  The records from `inp` are extracted
    and returned as a dictionary (that maps identifiers to (data,begin)
    pair).  If `mode` is 'anom' then data are converted to monthly
    anomalies; if `mode` is 'annual' then data are converted to annual
    anomalies (using the GISTEMP algorithm that copes with missing
    months).
    """

    # Clear Climate Code, tool directory
    import v2index
    # Clear Climate Code
    from code import series

    v2 = v2index.File(inp)

    table = {}
    for id in arg:
        for id12,rows in v2.get(id):
            data,begin = from_lines(rows, scale)
            if mode == 'anom':
                series.anomalize(data, None)
            if mode == 'annual':
                _, data = series.monthly_annual(data)
            table[id12] = (data,begin)

    return table

def colour_iter(background=(255,255,255)):
    """Generate a random sequence of colours, all different."""

    import random

    def distance(u, v):
        return sum((p-q)**2 for p,q in zip(u,v))**0.5

    # Internally all colour components are from 0 to 1; so convert
    # background.
    background = [c/255.0 for c in background]
    
    # Randomly generate colours and check that they are distance *b*
    # from the background colour and *r* from each other.  If we take
    # too many tries to find one... Make *r* smaller.

    b = 0.3
    r = 0.4
    # Ratio to reduce *r* by
    s = 0.7
    # Number of failures before reducing *r*.
    n = 8
    # List of all colours generated.
    l = []
    # Number of fails.
    fail = 0
    while 1:
        c = [random.random() for _ in range(3)]
        if (distance(c, background) > b and
          all(distance(c, x) > r for x in l)):
            yield "#%02x%02x%02x" % tuple(round(255*x) for x in c)
        fail += 1
        if fail > n:
            fail = 0
            r *= s


def update_config(config, v):
    """*config* is a configuration object used to store parameters.  *v*
    is an argument string of the form "parm1=value1;parm2=value2;...".
    Each "parm=value" pair sets an attribute of the config object.
    """

    l = v.split(';')
    for binding in l:
        attr,value = binding.split('=')
        attr = attr.strip()
        value = value.strip()
        for convert in [int, float, str]:
            try:
                value = convert(value)
                break
            except ValueError:
                pass
        setattr(config, attr, value)
    return config

def parse_topt(v):
    """Parse the t option which restricts the years to a particular
    range.  *v* is a string that is 2 (4-digit) years separated by a
    comma.  A pair of years is returned.
    """

    return map(int, v.split(','))

class Usage(Exception):
    pass
        
def main(argv=None):
    import codecs
    # http://www.python.org/doc/2.4.4/lib/module-getopt.html
    import getopt
    import sys
    if argv is None:
        argv = sys.argv

    try:
        infile = 'input/v2.mean'
        metafile = 'input/v2.inv'
        opt,arg = getopt.getopt(argv[1:], 'ac:o:d:m:s:t:y')
        if not arg:
            raise Usage('At least one identifier must be supplied.')
        outfile = arg[0] + '.svg'
        key = {}
        for k,v in opt:
            if k == '-c':
                update_config(config, v)
            if k == '-a':
                key['mode'] = 'anom'
            if k == '-o':
                outfile = v
            if k == '-d':
                infile = v
            if k == '-m':
                metafile = v
            if k == '-t':
                key['timewindow'] = parse_topt(v)
            if k == '-y':
                key['mode'] = 'annual'
            if k == '-s':
                key['scale'] = float(v)
        if outfile == '-':
            outfile = sys.stdout
        else:
            outfile = open(outfile, 'w')
        # See http://drj11.wordpress.com/2007/05/14/python-how-is-sysstdoutencoding-chosen/#comment-3770
        outfile = codecs.getwriter('utf-8')(outfile)
        if infile == '-':
            infile = sys.stdin
        else:
            infile = open(infile)
        metafile = open(metafile)
        derive_config(config)
        return plot(arg, inp=infile, out=outfile, meta=metafile, **key)
    except (getopt.GetoptError, Usage), e:
        sys.stdout.write('%s\n' % str(e))
        sys.stdout.write(__doc__)
        return 99

if __name__ == '__main__':
    sys.exit(main())

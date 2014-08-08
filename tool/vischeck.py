#!/usr/bin/env python
#
# vischeck.py - Visual Check
#
# David Jones, Climate Code Foundation
#
# Script to produce a visual check of the results (of STEP 5).

"""
vischeck.py [-o offset] GLB.txt [...]

Visually check a file of global anomalies, by converting it to a Google
Chart.  A URL is output (on stdout), visit/display that URL to see the
graph.

The file arguments should be in the GISTEMP tabular format, which is the
format used for this GISTEMP file
http://data.giss.nasa.gov/gistemp/tabledata/GLB.Ts.txt
and others.  The arguments can either be files on the local disk, or
URLs (in which case the data will be downloaded and used).

Multiple series (one per file) can be displayed on the same chart,
simply specify all the files on the command line (although more than 2
may run into problems with long URLs and the Google Chart API).

Two trend lines are drawn for each series: one for the entire series
and one for the last thirty years of data.

If series are very close, they will be displayed on top of each other.
An offset can be introduced between each series to shift it up the
chart.  -o offset will shift each series up (for a positive offset) or
down (negative offset); its units are the same as in the data files
(centikelvin usually).
"""

import datetime
import itertools
import math

class Error(Exception):
    """Some sort of problem."""

def annual_anomalies(f, extract=(65,72)):
    """
    Convert the text file *f* into a sequence of annual anomalies.
    An input file is expected to be one of the NH.*, SH.*, or GLB.*
    files that are the results of GISTEMP step 5.  The return value
    is an iterable over a sequence of pairs (year, datum); *datum*
    is an integer.  Years with invalid (missing) data are absent
    from the result stream.  Normally the data can be interpreted
    as centi-Kelvin.
    
    *extract* allows a different range of columns to be extracted for
    the anomaly data. It can either be a pair of numbers to
    specify the start and stop of a range (Python style, starts
    from 0, stop is not included), or it can be a string which
    is matched against the header line that starts "Year". Thus
    extract="JJA" will extract the series for the Northern
    summer season.
    """

    import re

    # Proceed by assuming that a line contains data if and only if it
    # starts with a 4-digit year; other lines are assumed to be
    # header/footer or decorative documentation.
    # This allows this function to work with both the direct output of
    # step 5, and also with the GISS published table data (which
    # includes more decorative documentation).
    for l in f:
        if l.startswith("Year") and isinstance(extract, basestring):
            # Use column numbers that are matched by *extract*.
            at = l.index(extract)
            end = l.index(extract) + len(extract)
            start = re.search(r' *$', l[:at]).start()
            extract = (start, end)
        if re.match(r'\d{4}', l):
            year = int(l[:4])
            try:
                yield (year, int(l[extract[0]:extract[1]]))
            except ValueError:
                pass

def trend(data):
    """Computes linear regression parameters (a,b) on the *data*.  The
    line y = a + b*x gives the best fit line to data.
    """
    sxx = sxy = syy = sx = sy = n = 0
    for (x,y) in data:
        if y is not None:
            sxx += x*x
            syy += y*y
            sx += x
            sy += y
            sxy += x * y
            n += 1
    if n < 2:
        return None,None,None
    # Make n a float. This contaminates all the subsequent divisions, making
    # them floating point divisions with floating point answers, which
    # is what we want.
    n = float(n)
    xbar = sx / n
    ybar = sy / n
    ssxx = sxx - (sx * sx) / n
    ssyy = syy - (sy * sy) / n
    ssxy = sxy - (sx * sy) / n
    b = ssxy / ssxx
    a = ybar - b * xbar
    r2 = (ssxy * ssxy) / (ssxx * ssyy)
    return (a,b, r2)

def asgooglechartURL(seq, options={}):
    """*seq* is a sequence of iterables (each one assumed to be the
    output from :meth:`asann`); convert into a URL suitable for passing
    to the Google Chart API.

    Each element of the sequence corresponds to a different data series,
    all the series are plotted on the same chart.

    *options* is a dict mapping keywords to optional attributes:

    options.offset: specifies the inter-chart offset, as per the -o option
    (see module docstring).  The default is zero.

    options.size: specifies the chart size as a tuple (width,
    height). The default is (600, 500).
    """

    # default options
    default_options = dict(offset=0,
                           size=(600,500),
                           colour=['ff0000', '000000', '0000ff', 'ff00ff'],
                           )
    default_options.update(options)
    options = default_options

    prefix = 'http://chart.apis.google.com/chart'

    # Read all data because we need to compute min and max years, and
    # min and max data values.
    data = [list(s) for s in seq]
    ymin,ymax = reasonable_yscale(data)

    yearmin = min([1880] + map(lambda p:p[0], itertools.chain(*data)))
    current_year = datetime.datetime.now().year
    yearmax = max([current_year] + map(lambda p:p[0], itertools.chain(*data)))
    data = [pad(s, yearmin, yearmax) for s in data]
    # Let y be the list of years for the chart legend.  We include
    # the first year of the series, the last year, and every decade
    # beginning.
    y = map(lambda x:x[0], data[0])
    for i in range(1,len(y)-1):
        if y[i]%10:
            y[i]=''
    xaxis = '|' + '|'.join(map(str, y))

    # Trendlines returns a trend object (describing the 2 trend lines
    # (short and long)); get one for each data series.
    trends = [trendlines(l) for l in data]
    trendfrags = [t.urlfrag for t in trends]

    # Arrange the series for the "chd=" parameter of the URL.
    # In an X-Y chart (cht=lxy) each series as drawn consists of a
    # sequence of xs and a sequence of ys.  In order we have...
    #   data series 1 (xs then ys),
    #   data series 2 (xs then ys),
    #   ...
    #   long trend 1 (xs then ys),
    #   short trend 1 (xs then ys),
    #   long trend 2 (xs then ys),
    #   short trend 2 (xs then ys)
    #   ...
    ds = ['-999|' + chartsingle(l) for l in data] + trendfrags

    vaxis = vaxis_labels(ymin,ymax)
    chxl = 'chxl=0:'+xaxis+'|1:'+vaxis+'|2:'+vaxis
    chxt = 'chxt=x,y,r'
    chd='chd=t:' + '|'.join(ds)
    chs='chs=' + 'x'.join(map(str, options['size']))
    # Choose scale, and deal with offset if we have to.
    # The scales are in the same order as the series (see above).  Each
    # series is associated with 4 scale values: xmin,xmax,ymin,ymax.
    offset = options['offset']
    scale = [-100,100,ymin,ymax] * (3 * len(seq))
    if offset and len(seq) > 1:
        for i in range(len(seq)):
            off = i*offset
            # Y range of actual data series
            scale[4*i+2] -= off
            scale[4*i+3] -= off
            # Y range of first trend
            scale[4*len(seq) + 8*i + 2] -= off
            scale[4*len(seq) + 8*i + 3] -= off
            # Y range of second trend
            scale[4*len(seq) + 8*i + 6] -= off
            scale[4*len(seq) + 8*i + 7] -= off
    chds = 'chds=' + ','.join(map(str, scale))
    # red, black, blue, magenta for the underlying graphs
    colours = options['colour']
    colours = colours[:len(data)]

    chm = 'chm=' + '|'.join(slope_markers(trends, colours))

    # Replicate each colour by three (for the chart and its two trend
    # lines).
    series_colours = colours + list(itertools.chain(*([c,c] for
      c in colours)))
    chco = 'chco=' + ','.join(series_colours)
    # Line Style
    chls = 'chls=' + '|'.join(['1']*len(data) + ['1,8,2', '1']*len(data))

    return prefix + '?' + '&'.join(
      ['cht=lxy',chds,chd,chxt,chxl,chco,chls,chm,chs])

def pad(data, yearmin, yearmax):
    """pad so that data series starts at yearmin and ends at yearmax."""

    nonelots = [None]*(yearmax-yearmin+1)

    if not data:
        return zip(range(yearmin, yearmax+1), nonelots)

    t0 = data[0][0]
    t1 = data[-1][0]

    assert yearmin <= t0 <= t1 <= yearmax
    return (zip(range(yearmin, t0), nonelots) +
      data +
      zip(range(t1+1, yearmax+1), nonelots))

def reasonable_yscale(data):
    """Examine the data and return a reasonable y scale as a min and max
    value."""

    tick = 50
    # Ensures subsequent divisions are in floating point.
    tick = float(tick)

    allvalues = [v for _,v in itertools.chain(*data) if v is not None]
    # Add values to the data to: ensure that a min and max can be
    # computed (otherwise empty list fails); and, to enforce a Y-axis
    # that covers at least -1 to +1.
    datamin = min(allvalues + [-100])
    datamax = max(allvalues + [+100])
    ymin = int(tick * math.floor(datamin/tick))
    ymax = int(tick * math.ceil(datamax/tick))
    return ymin,ymax

def vaxis_labels(ymin,ymax):
    """Return a reasonable vaxis label string for the google chart.
    Assumes that the data is in centikelvin."""

    # Spacing of ticks, native units.
    tick = 50

    assert ymin % tick == 0
    assert ymax % tick == 0
    l = ["%+.1f" % (v*0.01) for v in range(ymin,ymax+1,tick)]
    # Set first and last label to be blank.
    l[0] = '';
    l[-1] = '';

    return '|'.join([''] + l)

def slope_markers(trends, colours):
    """Create the markers to denote slopes / trends.  *trends* is a list
    of trend objects as returned by *trendlines*.  Returns a list."""

    import itertools
    import urllib

    l = [u'Trend in \N{DEGREE SIGN}C/century (R\N{SUPERSCRIPT TWO})'.encode('utf-8')]
    rowcolours = ['000000']
    for t,c in zip(trends,colours):
        if t.b_full is None:
            continue
        s = format_slope(['full','30-year'], [t.b_full, t.b_short],
          [t.r2_full, t.r2_short])
        l.append(s)
        rowcolours.append(c)
    return [
      urllib.quote_plus('@t%s,%s,0,%.2f:%.2f,12' % (
        text,
        colour,
        0.4, 0.2 - 0.05*row)) for text, colour, row in
        zip(l, rowcolours, itertools.count())]

def format_slope(texts, slopes, coefficients):
    """Return a string for the slopes and coefficients."""

    return '; '.join('%s: %.2f (%.2f)' % p for p in zip(texts, slopes, coefficients))

class Struct:
    pass

def trendlines(data):
    """Return a a triple of (url,slopelong,slopeshort) for
    the full and 30-year trend lines (url is a fragment).
    """

    result = Struct()
    # full trend
    (a,b,r2) = trend(data)
    if a is not None:
        yearmin = data[0][0]
        yearmax = data[-1][0]
        full_left_y = int(round(a + yearmin * b))
        full_right_y = int(round(a + yearmax * b))
        # thirty-year trend
        # Find most recent 30 years of _valid_ data
        valid_count = 0
        last_valid = -9999 # largest index with valid data
        for i,(x,y) in reversed(list(enumerate(data))):
            if y is None:
                continue
            last_valid = max(last_valid, i)
            valid_count += 1
            if valid_count >= 30:
                break
        (a_30,b_30,r2_30) = trend(data[i:])
        left_y = int(round(a_30 + (yearmin+i) * b_30))
        left_x = -100 + 200*float(i)/(yearmax-yearmin)
        right_y = int(round(a_30 + (yearmin+last_valid) * b_30))
        right_x = -100 + 200*float(last_valid)/(yearmax-yearmin)
        result.urlfrag = ("-100,100|%d,%d|%.0f,%.0f|%d,%d" %
            (full_left_y, full_right_y,
            left_x, right_x,
            left_y, right_y))
        result.b_full = b
        result.b_short = b_30
        result.r2_full = r2
        result.r2_short = r2_30
    else:
        result.b_full = None
        result.b_short = None
        result.r2_full = None
        result.r2_short = None
        result.urlfrag = "-999|-999|-999|-999"
    return result

def chartsingle(l):
    """Take a list and return a URL fragment for its Google
    chart."""

    d = map(lambda x:x[1], l)

    # Google Chart API says "Values less than the specified minimum are
    # considered missing values".  So we don't actually do any scaling
    # (we used to, but now it's all in the chds parameter).
    def scale(x):
        if x is None:
            return -999
        return x

    d = map(scale, d)
    return ','.join(map(str, d))

import sys

def chartit(fs, options={}, out=sys.stdout):
    """Convert the list of files *fs* to a Google Chart API url and print it
    on *out*.
    
    For documentation on *options* see `asgooglechartURL`.
    """

    import re
    import urllib

    k = {}
    if 'extract' in options:
        k['extract'] = options['extract']

    def anom(f):
        """Extract anomalies from file."""
        return annual_anomalies(f, **k)

    url = asgooglechartURL(map(anom, fs), options)
    if 'download' in options:
        img = urllib.urlopen(*url.split('?'))
        out.write(img.read())
    else:
        print >>out, url

def main(argv=None):
    import getopt
    # http://www.python.org/doc/2.4.4/lib/module-urllib.html
    import urllib

    if argv is None:
        argv = sys.argv

    options={}
    try:
        opt,arg = getopt.getopt(argv[1:], 'x:o:',
          ['extract=', 'offset=', 'size=', 'colour=', 'download'])
    except getopt.GetoptError, e:
        print >> sys.stderr, e.msg
        print >> sys.stderr, __doc__
        return 2
    for o,v in opt:
        if o in ('-o', '--offset'):
            options['offset'] = float(v)
        if o == '--size':
            options['size'] = v.split(',')
            if len(options['size']) != 2:
                raise Error("--size w,h is required")
        if o == '--colour':
            options['colour'] = v.split(',')
        if o in ('-x', '--extract'):
            try:
                options['extract'] = map(int, v.split(','))
            except ValueError:
                options['extract'] = v
        if o == '--download':
            options['download'] = True
            
    if not arg:
        print >> sys.stderr, __doc__
        return 2

    fs = map(urllib.urlopen, arg)
    chartit(fs, options)
    return 0

if __name__ == '__main__':
    sys.exit(main())

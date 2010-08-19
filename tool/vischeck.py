#!/usr/bin/env python
# $URL$
# $Rev$
#
# vischeck.py - Visual Check
#
# David Jones, Ravenbrook Limited
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

class Error(Exception):
    """Some sort of problem."""

def anomalies(f, extract=(65,72)):
    """Convert the text file *f* into a sequence of global anomalies.  An
    input file is expected to be one of the NH.*, SH.*, or GLB.* files
    that are the results of GISTEMP step 5.  The return value is an
    iterable over a sequence of pairs (year, datum); when present the
    datum is an integer, when not present it appears as None.  Normally
    the data can be interpreted as centi-Kelvin.  *extract* allows a
    different range of columns to be extracted for the anomaly data.
    """

    import re

    # Proceed by assuming that a line contains data if and only if it
    # starts with a 4-digit year; other lines are assumed to be
    # header/footer or decorative documentation.
    # This allows this function to work with both the direct output of
    # step 5, and also with the GISS published table data (which
    # includes more decorative documentation).
    for l in f:
        if re.match(r'\d{4}', l):
            year = int(l[:4])
            try:
                yield (year, int(l[extract[0]:extract[1]]))
            except ValueError:
                yield (year, None)

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
    output from :meth:`asann`) into a URL suitable for passing to the
    Google Chart API.

    Each element of the sequence corresponds to a different data series,
    all the series are plotted on the same chart.

    *Options* is a dict mapping keywords to optional attributes:

    option.offset: specifies the inter-chart offset, as per the -o option
    (see module docstring).  The default is zero.

    option.size: specifies the chart size as a tuple (width,
    height). The default is (600, 500).
    """

    import itertools

    # default options
    default_options = dict(offset=0,
                           size=(600,500),
                           colour=['ff0000', '000000', '0000ff', 'ff00ff'],
                           )
    default_options.update(options)
    options = default_options

    prefix = 'http://chart.apis.google.com/chart'

    # Read all data because we need to compute min and max years
    data = [list(s) for s in seq]

    yearmin = min(map(lambda p:p[0], itertools.chain(*data)))
    yearmax = max(map(lambda p:p[0], itertools.chain(*data)))
    data = [pad(s, yearmin, yearmax) for s in data]
    # Let y be the list of years for the chart legend.  We include
    # the first year of the series, the last year, and every decade
    # beginning.
    y = map(lambda x:x[0], data[0])
    for i in range(1,len(y)-1):
        if y[i]%10:
            y[i]=''
    xaxis = '|' + '|'.join(map(str, y))

    trends = [trendlines(l) for l in data]
    trenddata = [t[0] for t in trends]
    slopes = [t[1:len(data)+1] for t in trends]
    coefficients = [t[len(data)+1:] for t in trends]
    ds = ['-999|' + chartsingle(l) for l in data] + trenddata

    vaxis = '||-0.5|+0.0|+0.5|'
    chxl = 'chxl=0:'+xaxis+'|1:'+vaxis+'|2:'+vaxis
    chxt = 'chxt=x,y,r'
    chd='chd=t:' + '|'.join(ds)
    chs='chs=' + 'x'.join(map(str, options['size']))
    # Choose scale, and deal with offset if we have to.
    offset = options['offset']
    scale = [-100,100]*6
    if offset and len(seq) > 1:
        scale *= len(seq)
        for i in range(len(seq)):
            # Y range of actual data series
            scale[4*i+2] -= i*offset
            scale[4*i+3] -= i*offset
            # Y range of first trend
            scale[4*len(seq) + 8*i + 2] -= offset
            scale[4*len(seq) + 8*i + 3] -= offset
            # Y range of second trend
            scale[4*len(seq) + 8*i + 6] -= offset
            scale[4*len(seq) + 8*i + 7] -= offset
    chds = 'chds=' + ','.join(map(str, scale))
    # red, black, blue, magenta for the underlying graphs
    colours = options['colour']
    colours = colours[:len(data)]

    chm = 'chm=' + '|'.join(slope_markers(slopes, coefficients, colours))

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

    t0 = data[0][0]
    t1 = data[-1][0]

    assert yearmin <= t0 <= t1 <= yearmax
    nonelots = [None]*(yearmax-yearmin+1)
    return (zip(range(yearmin, t0), nonelots) +
      data +
      zip(range(t1+1, yearmax+1), nonelots))

def slope_markers(slopes, coefficients, colours):
    """Create the markers to denote slopes / trends."""

    import itertools
    import urllib

    l = [u'Trend (\N{DEGREE SIGN}C/Century) and R\N{SUPERSCRIPT TWO}'.encode('utf-8')] + [format_slope(['full','30-year'], *pl) for pl in zip(slopes, coefficients)]
    colours = ['000000'] + colours
    return [
      urllib.quote_plus('@t%s,%s,0,%.2f:%.2f,12' % (
        text,
        colour,
        0.4, 0.2 - 0.05*row)) for text, colour, row in
        zip(l, colours, itertools.count())]

def format_slope(texts, slopes, coefficients):
    """Return a string for the slopes and coefficients in pl, which is a list of pairs."""
    return ' / '.join('%s: %.2f (%.2f)' % p for p in zip(texts, slopes, coefficients))

def trendlines(data):
    """Return a a triple of (url,slopelong,slopeshort) for
    the full and 30-year trend lines (url is a fragment)."""
    # full trend
    (a,b,r2) = trend(data)
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
    return ("-100,100|%d,%d|%.0f,%.0f|%d,%d" % (full_left_y, full_right_y,
                                                left_x, right_x,
                                                left_y, right_y),
            b, b_30, r2, r2_30)

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
        return anomalies(f, **k)

    url = asgooglechartURL(map(anom, fs), options)
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
          ['offset=', 'size=', 'colour='])
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
        if o == '-x':
            options['extract'] = map(int, v.split(','))
            
    if len(arg):
        fs = map(urllib.urlopen, arg)
    else:
        fs = [sys.stdin]
    chartit(fs, options)
    return 0

if __name__ == '__main__':
    sys.exit(main())

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
Currently the same beginning and end years have to be used for all the
series.

If series are very close, they will be displayed on top of each other.
An offset can be introduced between each series to shift it up the
chart.  -o offset will shift each series up (for a positive offset) or
down (negative offset); its units are the same as in the data files
(centikelvin usually).
"""

class Error(Exception):
    """Some sort of problem."""

def asann(f):
    """Convert the text file *f* into a sequence of global anomalies.  An
    input file is expected to be one of the NH.*, SH.*, or GLB.* files
    that are the results of GISTEMP step 5.  The return value is an
    iterable over a sequence of pairs (year, datum); when present the
    datum is an integer, when not present it appears as None.  Normally
    the data can be interpreted as centi-Kelvin.
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
                yield((year, int(l[65:72])))
            except ValueError:
                yield((year, None))

def asgooglechartURL(seq, offset):
    """*seq* is a sequence of iterables (each one assumed to be the
    output from :meth:`asann`) into a URL suitable for passing to the
    Google Chart API.

    Each element of the sequence corresponds to a different data series,
    all the series are plotted on the same chart.

    *offset* corresponds to the -o option
    """

    prefix = 'http://chart.apis.google.com/chart'

    ds = []
    oy = None
    for i in seq:
        l = list(i)
        d = chartsingle(l)
        ds.append(d)
        # Let y be the list of years for the chart legend.  We include
        # the first year of the series, the last year, and every decade
        # beginning.
        y = map(lambda x:x[0], l)
        for i in range(1,len(y)-1):
            if y[i]%10:
                y[i]=''
        if oy is not None and oy != y:
            raise Error("Year ranges in data series are different.")
        oy = y

    xaxis = '|' + '|'.join(map(str, y))
    vaxis = '||-0.5|+0.0|+0.5|'
    chxl = 'chxl=0:'+xaxis+'|1:'+vaxis+'|2:'+vaxis
    chxt = 'chxt=x,y,r'
    chd='chd=t:' + '|'.join(ds)
    chs='chs=600x500'
    # Choose scale, and deal with offset if we have to.
    scale = [-100,100]
    if offset and len(seq) > 1:
        scale *= len(seq)
        for i in range(len(seq)):
            scale[2*i] -= i*offset
            scale[2*i+1] -= i*offset
    chds = 'chds=' + ','.join(map(str, scale))
    colours = ['ff0000','000000']
    chco = 'chco=' + ','.join(colours)
    return prefix + '?' + '&'.join(['cht=lc',chds,chd,chxt,chxl,chco,chs])

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

def chartit(fs, offset):
    """Convert the list of files *fs* to a Google Chart API url and print it
    on standard output.
    
    If the list is exactly one element long then a PNG file received
    from the URL is written to a file 'foo.png' where its name is
    derived from the input file.

    *offset* specifies the inter-chart offset, as per the -o option (see
    module docstring).
    """

    import re
    import urllib

    url = asgooglechartURL(map(asann, fs), offset=offset)
    print url
    uf = urllib.urlopen(url)
    if len(fs) != 1:
        return
    f = fs[0]

    def name(f):
        """Return the "name" of a file object.  Works with genuine
        objects and objects return from urllib.urlopen.
        """

        try:
            return f.name
        except:
            return f.geturl()
    # Open a PNG output file with a name based on the input.
    # 'foo.txt' -> 'foo.png' ; '<stdin>' -> '<stdin>.png'
    o = open(re.sub('(\.[^.]*|)$', '.png', name(f)), 'wb')
    while True:
        x = uf.read(1000)
        if x == '':
            break
        o.write(x)
    o.close()

import sys

def main(argv=None):
    import getopt
    # http://www.python.org/doc/2.4.4/lib/module-urllib.html
    import urllib

    if argv is None:
        argv = sys.argv

    # offset between each series (-o option)
    offset = 0

    try:
        opt,arg = getopt.getopt(argv[1:], 'o:', ['offset='])
    except getopt.GetoptError, e:
        print >> sys.stderr, e.msg
        print >> sys.stderr, __doc__
        return 2
    for o,v in opt:
        if o in ('-o', '--offset'):
            offset = float(v)
    if len(arg):
        fs = map(urllib.urlopen, arg)
    else:
        fs = [sys.stdin]
    chartit(fs, offset=offset)
    return 0

if __name__ == '__main__':
    sys.exit(main())

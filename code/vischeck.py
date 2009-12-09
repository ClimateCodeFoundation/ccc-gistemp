#!/usr/bin/env python
# $URL$
# $Rev$
#
# vischeck.py - Visual Check
#
# David Jones, Ravenbrook Limited
#
# Script to produce a visual check of the results (of STEP 5).

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

def asgooglechartURL(seq):
    """*seq* is a sequence of iterables (each one assumed to be the
    output from :meth:`asann`) into a URL suitable for passing to the
    Google Chart API.

    Each element of the sequence corresponds to a different data series,
    all the series are plotted on the same chart.
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
    chds = 'chds=-100,100'
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

def chartit(fs):
    """Convert the list of files *fs* to a Google Chart API url and print it
    on standard output.
    
    If the list is exactly one element long then a PNG file received
    from the URL is written to a file 'foo.png' where its name is
    derived from the input file.
    """

    import re
    import urllib

    url = asgooglechartURL(map(asann, fs))
    print url
    uf = urllib.urlopen(url)
    if len(fs) != 1:
        return
    f = fs[0]
    # Open a PNG output file with a name based on the input.
    # 'foo.txt' -> 'foo.png' ; '<stdin>' -> '<stdin>.png'
    o = open(re.sub('(\.[^.]*|)$', '.png', f.name), 'wb')
    while True:
        x = uf.read(1000)
        if x == '':
            break
        o.write(x)
    o.close()

def main(argv=None):
    import sys

    if argv is None:
        argv = sys.argv
    if len(argv[1:]):
        fs = map(open, argv[1:])
    else:
        fs = [sys.stdin]
    chartit(fs)

if __name__ == '__main__':
    main()

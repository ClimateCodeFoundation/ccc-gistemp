#!/usr/bin/env python
# $URL$
# $Rev$
#
# vischeck.py - Visual Check
#
# David Jones, Ravenbrook Limited
#
# Script to produce a visual check of the results (of STEP 5).

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

def asgooglechartURL(i):
    """Convert an iterable (assumed to be output from :meth:`asann`) into
    a URL suitable for passing to the Google Chart API.
    """

    prefix = 'http://chart.apis.google.com/chart'

    l = list(i)

    d = map(lambda x:x[1], l)

    # Google Chart API says you can specify an out-of-range value for
    # a missing datum, but this doesn't work for line charts if data
    # scaling is being used.  So we drop the data scaling, scaling it
    # ourselves to a range of 0 to 100, and use -1 for the missing value.
    def scale(x):
        if x is None:
            return -1
        x = (x + 100.0) / 2.0
        if not (0 <= x <= 100):
            x = -1
        return x

    d = map(scale, d)
    y = map(lambda x:x[0], l)
    for i in range(1,len(y)-1):
        if y[i]%10:
            y[i]=''
    xaxis = '|' + '|'.join(map(str, y))
    vaxis = '||-0.5|+0.0|+0.5|'
    chxl = 'chxl=0:'+xaxis+'|1:'+vaxis+'|2:'+vaxis
    chxt = 'chxt=x,y,r'
    chd='chd=t:' + ','.join(map(str, d))
    chs='chs=600x500'
    chco='chco=ff0000'
    return prefix + '?' + '&'.join(['cht=lc',chd,chxt,chxl,chco,chs])

def chartit(f):
    """Convert the file object *f* to a Google Chart API url and print it
    on standard output, and *ahem* write the PNG received from the same
    URL to a file 'foo.png' where its name is derived from the input
    file.
    """

    import re
    import urllib

    url = asgooglechartURL(asann(f))
    print url
    uf = urllib.urlopen(url)
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
    for f in fs:
        chartit(f)

if __name__ == '__main__':
    main()

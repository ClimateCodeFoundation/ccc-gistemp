#!/usr/bin/env python
# $URL$
# $Id$
# 
# compare_results.py -- compare results from two runs of the GISTEMP algorithm
#
# Gareth Rees, 2009-12-17

"""compare_results.py [options] RESULTA RESULTB

Compare results from two runs of the GISTEMP algorithm.  Specify two
directories RESULTA and RESULTB on the command line: each should be a
directory containing the results of a run of the GISTEMP algorithm.  A
comparison of the two sets of results is output to the file index.html
in the current directory.

Options:
   --help          Print this text.
   --output=FILE   Write the output to this file (default: index.html).
   --labela=LABEL  Label for RESULTA (default: CCC).
   --labelb=LABEL  Label for RESULTB (default: GISTEMP).
"""

# http://www.python.org/doc/2.4.4/lib/module-getopt.html
import getopt
# http://www.python.org/doc/2.4.4/lib/module-math.html
import math
# http://www.python.org/doc/2.4.4/lib/module-os.html
import os
# http://www.python.org/doc/2.4.4/lib/module-re.html
import re
# http://www.python.org/doc/2.4.4/lib/module-sys.html
import sys
# http://www.python.org/doc/2.4.4/lib/module-xml.sax.saxutils.html
import xml.sax.saxutils
# Clear Climate Code
import vischeck

class Fatal(Exception):
    def __init__(self, msg):
        self.msg = msg

# This is derived from the similar function vischeck.asann.
def asmon(f):
    """Convert the text file *f* into a sequence of monthly anomalies.
    An input file is expected to be one of the NH.*, SH.*, or GLB.*
    files that are the results of GISTEMP step 5.  The return value is
    an iterable over a sequence of tuples (year, month, datum) with
    month running from 0 to 11; when present the datum is an integer,
    when not present it appears as None.  Normally the data can be
    interpreted as centi-Kelvin.
    """

    # Proceed by assuming that a line contains data if and only if it
    # starts with a 4-digit year; other lines are assumed to be
    # header/footer or decorative documentation.
    # This allows this function to work with both the direct output of
    # step 5, and also with the GISS published table data (which
    # includes more decorative documentation).
    for l in f:
        if re.match(r'\d{4}', l):
            year = int(l[:4])
            for month in range(12):
                try:
                    yield((year, month, int(l[(month + 1) * 5:(month + 2) * 5])))
                except ValueError:
                    yield((year, month, None))

def stats(seq):
    """Return the mean and standard deviation of the numbers in the
    non-empty sequence *seq*."""
    assert seq
    n = float(len(seq))
    mean = sum(seq) / n
    variance = sum(map(lambda x: (x - mean) ** 2, seq)) / n
    return mean, math.sqrt(variance)

def difference(seqs, k, scale):
    """Return a sequence of differences between the *k*th items of the
    tuples in the two sequences *seqs*, multiplied by *scale*."""
    assert len(seqs) == 2
    def diff(a, b):
        if a[k] is None or b[k] is None:
            return 0
        else:
            return (a[k] - b[k]) * scale
    return map(diff, *seqs)

def distribution_url(d, w):
    """Return the URL of a Google chart showing the distribution of the
    values in the sequence *d*, collected into bins of width *w*."""

    # Map from integer multiple of w to number of values in that bin.
    bins = {}
    for v in d:
        b = int(math.floor(0.5 + v / w))
        bins[b] = bins.get(b, 0) + 1
    bin_min = min(bins.keys())
    bin_max = max(bins.keys())

    # Number of values in each bin.
    d = map(lambda x: bins.get(x, 0), range(bin_min, bin_max + 1))
    dmax = max(d)

    # Google chart documentation at http://code.google.com/apis/chart/
    rx = '%f,%f' % (bin_min * w, bin_max * w)
    ry = '0,%d,%d' % (dmax, 10 ** math.floor(math.log10(dmax)))
    chart = [
        'chs=600x300', # Chart size (pixels).
        'cht=bvg',     # Chart type: Horizontal bar chart, stacked bars.
        'chd=t:' + ','.join(map(str, d)),     # Data.
        'chds=0,%d' % dmax,                   # Data scaling. 
        'chxt=x,y,r',                         # Axes to label.
        'chxr=0,%s|1,%s|2,%s' % (rx, ry, ry), # Axes ranges.
        ]

    return 'http://chart.apis.google.com/chart?' + '&'.join(chart)

def differences_url(anns):
    """Return the URL of a Google chart showing the difference between
    the two sequences of annual anomalies in the sequence *anns*."""
    assert len(anns) == 2

    # Take difference.
    d = difference(anns, 1, 0.01)

    # List of years
    years = map(lambda a: a[0], anns[0])
    assert years == map(lambda a: a[0], anns[1])

    # Google chart documentation at http://code.google.com/apis/chart/
    dmin = min(d)
    dmax = max(d)
    ry = '%.2f,%.2f,.01' % (dmin, dmax)
    rx = '%d,%d,10' % (min(years), max(years))
    chart = [
        'chs=600x200',   # Chart size (pixels).
        'cht=bvg',       # Chart type: Horizontal bar chart, stacked bars.
        'chbh=a',        # Automatically resize bars to fit.
        'chd=t:%s' % ','.join(map(str, d)),   # Data values.
        'chds=%.2f,%.2f' % (dmin, dmax),      # Data scaling. 
        'chxt=x,y,r',                         # Axes to label.
        'chxr=0,%s|1,%s|2,%s' % (rx, ry, ry), # Axes ranges.
        ]

    return 'http://chart.apis.google.com/chart?' + '&'.join(chart)

def compare(dirs, labels, o):
    """Compare results from the two directories in the sequence *dirs*,
    which are labelled by the two strings in the sequence *labels*. Write a
    comparison document in HTML to the stream *o*."""
    assert len(dirs) == 2
    assert len(labels) == 2

    # XML-encode meta-characters &,<,>
    escape = xml.sax.saxutils.escape
    labels = map(escape, labels)

    title = "Comparison of %s and %s" % tuple(dirs)
    print >>o, """<!doctype HTML>
<html>
<head>
<title>%s</title>
</head>
<body>
<h1>%s</h1>
""" % (title, title)

    anomaly_file = '%s.Ts.ho2.GHCN.CL.PA.txt'
    for region, code in [
        ('Global',              'GLB'), 
        ('Northern hemisphere', 'NH'),
        ('Southern hemisphere', 'SH'),
    ]:
        # Annual series
        fs = map(lambda d: open(os.path.join(d, anomaly_file % code), 'r'), dirs)
        anns = map(list, map(vischeck.asann, fs))
        url = vischeck.asgooglechartURL(anns, offset=0)
        print >>o, '<h2>%s annual temperature anomaly</h2>' % region
        print >>o, '<p>%s in red, %s in black.</p>' % tuple(labels) 
        print >>o, '<img src="%s">' % escape(url)
        print >>o, ('<h3>%s annual residues (%s - %s)</h3>'
                    % (region, labels[0], labels[1]))
        print >>o, '<img src="%s">' % escape(differences_url(anns))

        diffs = difference(anns, 1, 0.01)
        print >>o, '<h3>%s annual residue distribution</h3>' % region
        print >>o, '<img src="%s">' % escape(distribution_url(diffs, 0.01))
        print >>o, '<h3>%s annual residue summary</h3>' % region
        print >>o, '<ul>'
        print >>o, '<li>Min = %f<li>Max = %f' % (min(diffs), max(diffs))
        print >>o, '<li>Mean = %f<li>Standard deviation = %f' % stats(diffs)
        print >>o, '</ul>'

        # Monthly series
        fs = map(lambda d: open(os.path.join(d, anomaly_file % code), 'r'), dirs)
        mons = map(asmon, fs)
        diffs = difference(mons, 2, 0.01)
        print >>o, '<h3>%s monthly residue distribution</h3>' % region
        print >>o, '<img src="%s">' % escape(distribution_url(diffs, 0.01))
        print >>o, '<h3>%s monthly residue summary</h3>' % region
        print >>o, '<ul>'
        print >>o, '<li>Min = %f<li>Max = %f' % (min(diffs), max(diffs))
        print >>o, '<li>Mean = %f<li>Standard deviation = %f' % stats(diffs)
        print >>o, '</ul>'

    print >>o, "</body>"
    print >>o, "</html>"
    o.close()

def main(argv = None):
    if argv is None:
        argv = sys.argv
    try:
        # Parse command-line arguments.
        output = 'index.html'
        labels = ['CCC', 'GISS']
        try:
            opts, args = getopt.getopt(argv[1:], 'ho:l:m:',
                                       ['help', 'output=', 'labela=', 'labelb='])
            for o, a in opts:
                if o in ('-h', '--help'):
                    print __doc__
                    return 0
                elif o in ('-o', '--output'):
                    output = a
                elif o in ('-l', '--labela'):
                    labels[0] = a
                elif o in ('-m', '--labelb'):
                    labels[1] = a
                else:
                    raise Fatal("Unsupported option: %s" % o)
            if len(args) != 2:
                raise Fatal("Expected two result directories, but got '%s'"
                            % ' '.join(args))
        except getopt.error, msg:
            raise Fatal(str(msg))

        # Do the comparison.
        compare(args, labels, open(output, 'w'))
        return 0
    except Fatal, err:
        sys.stderr.write(err.msg)
        sys.stderr.write('\n')
        return 2

if __name__ == '__main__':
    sys.exit(main())

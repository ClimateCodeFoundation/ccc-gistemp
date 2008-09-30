#!/usr/bin/env python
# $Id: //info.ravenbrook.com/project/ccc/master/code/subboxcmp.py#3 $
#
# subboxcmp.py
#
# Compare two binary gridded subbox files (as output from STEP3), often called
# SBBX1880.Ts.GHCN.CL.PA.1200
#
# David Jones.  Ravenbrook Limited.
#
# See subboxtotext.py for notes on the file format.
#
# The two files are assumed to be created from aligned (possibly the
# same) sources.  They are expected to contain the same number of
# subboxes, with the subboxes in the order and having the same
# boundaries.  Both files are expected to have the same: start year,
# averaging method (mavg field), BAD flag.  If any of these expectations
# is not met the comparison will be terminated immediately.
#
# The following differences are reported:
# differing number of stations combined for a given subbox;
# differing number of station months used for a given subbox;
# "large" differences in pointwise comparison of time series (the
#    difference that is reported is configurable and currently defaults
#    to 1e-4)
# "large" differences in d, the distance to nearest station.  Where
#    large means more than 0.5
# Record numbers start at 0 for the first record in the file which is a
# header.  Thus the records corresponding to subboxes are from 1 to 8000

import fort
import struct
import sys

def cmp(a, b, dt=1e-4, dd=0.5, output=sys.stdout, error=sys.stderr):
    """Compare two files."""

    # Compute the width of a standard word according to Python's struct
    # module...
    w = len(struct.pack('=I', 0))
    # and a suitable string format.
    # http://www.python.org/doc/2.3.5/lib/typesseq-strings.html
    # The string format is of the form '%08x' but the value of 8 may be
    # replaced.
    fmt = '%%0%dx' % (2*w)
    # Width of a float
    wf = len(struct.pack('f', 0.0))

    a = fort.File(a)
    b = fort.File(b)

    ra = a.readline()
    rb = b.readline()
    rn = 0

    # Number of words in header, preceding title.
    n = 8
    if ra[:n*w] != rb[:n*w]:
        error.write('headers differ:\n' +
            str(struct.unpack('%di' % n, ra[:n*w])) +
            str(struct.unpack('%di' % n, rb[:n*w])) +
            '\n')
        sys.exit(4)

    dmax = -1
    dmaxrn = None
    tmax = -1
    tmaxrni = ()

    while True:
        ra = a.readline()
        rb = b.readline()
        rn += 1
        if ra == None and rb == None:
            break
        if ra == None or rb == None:
            error.write('files differ in size')
            sys.exit(4)
        if len(ra) != len(rb):
            error.write('Record %d is different size' % rn)
            break

        traila = ra[-7*w:]
        trailb = rb[-7*w:]
        if traila[:4*w] != trailb[:4*w]:
            error.write('Record %d is for different boxes' % rn)
            break
        counta = struct.unpack('2I', traila[4*w:6*w])
        countb = struct.unpack('2I', trailb[4*w:6*w])
        if counta[0] != countb[0]:
            output.write('Record %d NSTNS: %d %d\n' %
                (rn, counta[0], countb[0]))
        if counta[1] != countb[1]:
            output.write('Record %d NSTMNS: %d %d\n' %
                (rn, counta[1], countb[1]))
        da = struct.unpack('f', traila[6*w:])[0]
        db = struct.unpack('f', trailb[6*w:])[0]
        if abs(da-db) >= dd:
            output.write('Record %d D: %s %s\n' %
                (rn, repr(da), repr(db)))
        if abs(da-db) >= dmax:
            dmax = abs(da-db)
            dmaxrn = rn
        
        ra = ra[:-7*w]
        rb = rb[:-7*w]
        n = len(ra)//wf  # number of time series entries
        ta = struct.unpack('%df' % n, ra)
        tb = struct.unpack('%df' % n, rb)
        for i in range(n):
            d = abs(ta[i]-tb[i])
            if d >= dt:
                output.write('Record %d data %i: %s %s diff: %s\n' %
                    (rn, i, repr(ta[i]), repr(tb[i]), repr(d)))
            if d >= tmax:
                tmax = d
                tmaxrni = (rn, i)
    output.write('Maximum difference in d (record %d): %s\n' %
        (dmaxrn, repr(dmax)))
    output.write('Maximum difference in t (record %d item %d): %s\n' %
        (tmaxrni[0], tmaxrni[1], repr(tmax)))

def main():
    cmp(open(sys.argv[1], 'rb'), open(sys.argv[2], 'rb'))

if __name__ == '__main__':
    main()

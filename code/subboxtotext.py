#!/usr/bin/env python
# $Id: //info.ravenbrook.com/project/ccc/master/code/subboxtotext.py#7 $
#
# subboxtotext.py
#
# Converts the gridded subbox file output from STEP3, often called
# SBBX1880.Ts.GHCN.CL.PA.1200, to a text format so that they can be
# compared.
#
# David Jones.  Ravenbrook Limited.
#
# The comments at the top of to.SBBXgrid.f that describe the output
# file formats are incorrect.  The actual format of this file is a
# Fortran binary file with an initial header record followed by a series
# of identically formatted records, one for each subbox:
#
# Record 1: Header
# Record 2: AVG array, lats, latn, lonw, lone, nstns, nstmns, d
# Record 3: ...
#
# The header comprises of 8 words followed by 80 characters:
# INFO(1)   1 (not trimmed)
#     (2)   KQ
#     (3)   MAVG
#     (4)   MONM (length of each time series)
#     (5)   reclen (length of each record in words)
#     (6)   IYRBEG (first year of time series)
#     (7)   BAD (bad data value, as integer)
#     (8)   precipitation trace flag
# TITLE     80-character string.
#
# Each subbox record consists of:
# MONM single-precision floats;
# 4 (signed) integers for the bounds of the subbox
# nstns is the number of stations combined into the subbox time series
# nstmn is the total number of data used from all the combined stations
# d is an approximation to the distance from the box centre of the
# nearest station used.

import fort
import struct
import sys

def totext(file, output=sys.stdout, error=sys.stderr, metaonly=False):
    """Convert binary gridded subbox file to text format.
    
    If metaonly is True then only the subbox metadata is output, the
    time series are not."""

    # Compute the width of a standard word according to Python's struct
    # module...
    w = len(struct.pack('=I', 0))
    # Width of a float
    wf = len(struct.pack('f', 0.0))
    # and a suitable string format.
    # http://www.python.org/doc/2.3.5/lib/typesseq-strings.html
    # The string format is of the form '%08x' but the value of 8 may be
    # replaced.
    fmt = '%%0%dx' % (2*w)

    f = fort.File(file)

    r = f.readline()

    # Number of words in header, preceding title.
    n = 8
    info = struct.unpack('%di' % n, r[:n*w])
    print info
    print r[n*w:]
    yrbeg = info[5]
    mavg = info[2]
    km = 1
    if mavg == 6:
        km = 12

    for r in f:
        trail = r[-7*w:]
        box = struct.unpack('4i', trail[:4*w])
        box = tuple(map(lambda x: x/100.0, box))
        nstns, nstmns = struct.unpack('2I', trail[4*w:6*w])
        d = struct.unpack('f', trail[6*w:])[0]
        loc = '%+06.2f%+06.2f%+07.2f%+07.2f' % box
        output.write('%s META %3d %6d %f\n' % (loc, nstns, nstmns, d))
        if metaonly:
            continue
        
        r = r[:-7*w]
        n = len(r)//wf  # number of time series entries
        t = struct.unpack('%df' % n, r)
        # 12 entries per output line, which is usually one year's worth,
        # (but see km).
        p = 12
        for i in range(len(t)//p) :
            output.write('%s %d' % (loc, (yrbeg + i*p//km)))
            for j in range(p) :
                k = i*p+j
                output.write(' %s' % repr(t[k]))
            output.write('\n')

def main(argv=None):
    import getopt

    if argv is None:
        argv = sys.argv

    metaonly = False
    opt, arg = getopt.getopt(argv[1:], 'm')
    for o,v in opt:
        if o == '-m':
            metaonly = True
    if len(arg) == 0:
        totext(sys.stdin, metaonly=metaonly)
    else:
        for n in arg:
            totext(open(n, 'rb'), metaonly=metaonly)

if __name__ == '__main__':
    main()

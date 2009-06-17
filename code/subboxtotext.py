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

def totext(file, output=sys.stdout, log=sys.stderr, metaonly=False,
  bos='>', trimmed='fromfile'):
    """Convert binary gridded subbox file to text format.
    
    If metaonly is True then only the subbox metadata is output, the
    time series are not.
    
    `trimmed` determines whether the input file is trimmed; it can be
    one of ``True`` (file is trimmed), ``False`` (file is not trimmed),
    or ``'fromfile'`` (will examine file to determine if it is trimmed
    or not).
    """

    assert trimmed in [True, False, 'fromfile']

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

    f = fort.File(file, bos=bos)

    r = f.readline()

    # Number of words in header, preceding title.
    n = 8
    info = struct.unpack(bos + ('%di' % n), r[:n*w])
    print info
    print r[n*w:]
    yrbeg = info[5]
    mavg = info[2]
    km = 1
    if mavg == 6:
        km = 12
    if trimmed == 'fromfile':
        trimmed = info[0] != 1
        print >>log, "Determined that trimmed=%s from file." % trimmed

    for r in f:
        if trimmed:
            meta = r[1*w:8*w]
            r = r[8*w:]
        else:
            meta = r[-7*w:]
            r = r[:-7*w]
        box = struct.unpack(bos + '4i', meta[:4*w])
        box = tuple(map(lambda x: x/100.0, box))
        nstns, nstmns = struct.unpack(bos + '2I', meta[4*w:6*w])
        d = struct.unpack(bos + 'f', meta[6*w:])[0]
        loc = '%+06.2f%+06.2f%+07.2f%+07.2f' % box
        n = len(r)//wf  # number of time series entries
        output.write('%s META %6d %3d %6d %f\n' % (loc, n, nstns, nstmns, d))
        if metaonly:
            continue
        
        t = struct.unpack(bos + ('%df' % n), r)
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
    bos = '>'
    trimmed = 'fromfile'
    opt, arg = getopt.getopt(argv[1:], 'b:mtu')
    for o,v in opt:
        if o == '-m':
            metaonly = True
        if o == '-b':
            bos = v
        if o == '-t':
            trimmed = True
        if o == '-u':
            trimmed = False
    if len(arg) == 0:
        totext(sys.stdin, metaonly=metaonly, bos=bos, trimmed=trimmed)
    else:
        for n in arg:
            totext(open(n, 'rb'), metaonly=metaonly, bos=bos,
              trimmed=trimmed)

if __name__ == '__main__':
    main()

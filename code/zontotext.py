#!/usr/bin/env python
# $URL$
# $Rev$
#
# zontotext.py
#
# Converts ZON.* file to plain text.
# 
# Nick Barnes.  Ravenbrook Limited.

import fort
import struct
import sys

# Not sure what these constants are; cribbed from zonav.f.
# They boil down to the fact that the ZON.* file contains 14 records.
# 

jbm = 8 # bands 90-64-44-24-0-24-44-64-90
nzs = 3 # zones 90-24-24-90
jzm = jbm + nzs + 3 # NH, SH, global

class Error(Exception):
    """Some problem."""

def totext(file, output=sys.stdout, log=sys.stderr, metaonly=False, bos='>'):
    """Convert zonal monthly averages to text format.
    
    If metaonly is True then only the zonal metadata is output, the
    time series are not.
    """

    # :todo: move into common module
    from zonav import swaw

    # The width of a standard word according to Python's struct module...
    w = len(struct.pack('=I', 0))

    f = fort.File(file, bos=bos)
    r = f.readline()

    # Number of words in header, preceding title.
    n = 8
    info = struct.unpack(bos + ('%di' % n), r[:n*w])
    output.write(repr(info))
    output.write('\n%s\n' % r[n*w:n*w+80])
    output.write('%s\n' % r[n*w+80:])

    # m: time frames per year
    if info[2] == 6:
        m = 12
    else:
        m = 4
    first_year = info[5]
    months = info[3]
    years = months/m
    last_year = first_year + years - 1

    # Each line contains N ar values and N weight (area?) values,
    # followed by an 80-character title string.
    descriptor = bos + '%df80s' % (months*2)

    for i in range(jzm):
        r = f.readline()
        if r is None:
            raise Error('Unexpected end of file.')
        data = struct.unpack(descriptor, r)
        output.write(swaw(data[-1]))
        if metaonly:
            continue
        for set in range(2):
            for year in range(first_year, last_year+1):
                offset = (year-first_year)*m + (months * set)
                output.write('%s[%4d]: %s\n' % (['AR','WT'][set],
                                                year,
                                                ' '.join(map(repr, data[offset:offset+m]))))

def main(argv=None):
    import getopt

    if argv is None:
        argv = sys.argv

    metaonly = False
    bos = '>'
    opt, arg = getopt.getopt(argv[1:], 'b:m')
    for o,v in opt:
        if o == '-m':
            metaonly = True
        if o == '-b':
            bos = v
    if len(arg) == 0:
        totext(sys.stdin, metaonly=metaonly, bos=bos)
    else:
        for n in arg:
            totext(open(n, 'rb'), metaonly=metaonly, bos=bos)

if __name__ == '__main__':
    main()

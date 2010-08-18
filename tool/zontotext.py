#!/usr/bin/env python
# $URL$
# $Rev$
#
# zontotext.py
#
# Converts ZON.* file to plain text.
# 
# Nick Barnes, David Jones.  Ravenbrook Limited.

import struct
import sys

# Clear Climate Code
import fort
import giss_io

# Not sure what these constants are; cribbed from zonav.f.
# They boil down to the fact that the ZON.* file contains 14 records.
# 

jbm = 8 # bands 90-64-44-24-0-24-44-64-90
nzs = 3 # zones 90-24-24-90
jzm = jbm + nzs + 3 # NH, SH, global

class Error(Exception):
    """Some problem."""

def totext(file, output=sys.stdout, log=sys.stderr, metaonly=False,
  bos='>', format='v2'):
    """Convert zonal monthly averages to text format.
    
    If metaonly is True then only the zonal metadata is output, the
    time series are not.
    """

    # The width of a standard word according to Python's struct module...
    w = len(struct.pack('=I', 0))

    f = fort.File(file, bos=bos)
    r = f.readline()

    # Number of words in header, preceding title.
    n = 8
    info = struct.unpack(bos + ('%di' % n), r[:n*w])
    if 'v2' != format:
        output.write(repr(info))
        output.write('\n%s\n' % r[n*w:n*w+80])
        output.write('%s\n' % r[n*w+80:])

    if 'v2' == format:
        v2out = giss_io.V2MeanWriter(file=output, scale=0.01)

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
        title = data[-1]
        if format == 'v2':
            title = v2title(title)
        else:
            output.write(title + '\n')
        if metaonly:
            continue
        for idx in range(2):
            if 'v2' == format and idx > 0:
                # Only output temps, not weights. :todo: fix this.
                continue
            for year in range(first_year, last_year+1):
                offset = (year-first_year)*m + (months * idx)
                temps = data[offset:offset+m]
                if 'v2' == format:
                    assert 12 == m
                    v2out.writeyear(title+('TW'[idx]), year, temps)
                else:
                    output.write('%s[%4d]: %s\n' %
                      (['AR','WT'][idx],
                      year,
                      ' '.join(map(repr, temps))))

def v2title(s):
    """Convert a title as it appears in the ZON file into an 11
    character V2 Mean style station identifier.

    The returned strings will be of the form "Z+SS.S+NN.N".
    """

    # http://docs.python.org/release/2.4.4/lib/module-re.html
    import re

    s = s.lower()

    # Three special cases for hemispheres and global...
    if 'global' in s:
        return "Z-90.0+90.0"
    if 'hemisphere' in s:
        if 'north' in s:
            return "Z+00.0+90.0"
        if 'south' in s:
            return "Z-90.0+00.0"

    # One of the "regular" cases. We extract the zone boundaries from
    # the title.

    s = s.replace('equator', ' 00.0 n')
    m = re.search(r'(\d\d\.\d\s+[ns]).*(\d\d\.\d\s+[ns])', s)
    if not m:
        raise ValueError("Cannot convert zone title string: %s" % s)
    bound = m.groups()
    def cvt(x):
        """Convert "DD.D n" to "+DD.D" and "DD.D s" to "-DD.D"."""

        x = x.replace('s', '-')
        x = x.replace('n', '+')
        return x[-1] + x[:4]
    return 'Z' + ''.join(cvt(b) for b in bound)

def main(argv=None):
    import getopt

    if argv is None:
        argv = sys.argv

    k = {}
    opt, arg = getopt.getopt(argv[1:], 'b:mt')
    for o,v in opt:
        if o == '-m':
            k['metaonly'] = True
        if o == '-b':
            k['bos'] = v
        if o == '-t':
            k['format'] = 'text'
            
    if len(arg) == 0:
        totext(sys.stdin, **k)
    else:
        for n in arg:
            totext(open(n, 'rb'), **k)

if __name__ == '__main__':
    main()

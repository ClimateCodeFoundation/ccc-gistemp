#!/usr/bin/env python
# $URL$
# $Rev$
#
# zontotext.py
# 
# Nick Barnes, David Jones.  Climate Code Foundation.

"""
Converts (fortran binary) ZON.* file to GHCN v2 format (which is a plain
text format).  Can also convert BX.* files.
"""


import struct
import sys

# Clear Climate Code
import extend_path
from code import eqarea
import fort
import gio

class Error(Exception):
    """Some problem."""

def totext(inp, output=sys.stdout, log=sys.stderr, metaonly=False,
  bos='>', format='v2'):
    """Convert monthly averages to text format; *inp* is the input file
    and should be either a binary zone file (ZON.*) or a binary box file
    (BX.*).
    
    If metaonly is True then only the zonal metadata is output, the
    time series are not.
    """

    # The width of a standard word according to Python's struct module...
    w = len(struct.pack('=I', 0))

    f = fort.File(inp, bos=bos)
    r = f.readline()

    # Number of words in header, preceding title.
    n = 8
    info = struct.unpack(bos + ('%di' % n), r[:n*w])
    title = r[n*w:]
    if 'zones' in title.lower():
        content = 'zones'
    else:
        content = 'boxes'

    if metaonly or 'v2' != format:
        output.write(repr(info))
        output.write('\n' + title + '\n')
    if metaonly:
        return
    if 'v2' == format:
        v2out = gio.GHCNV2Writer(file=output, scale=0.01)

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
    # followed by...
    #  - (for zones) an 80-character title string.
    #  - (for boxes) 5 words.
    rest = dict(zones=80, boxes=5*w)[content]
    descriptor = bos + '%df%ds' % (months*2, rest)

    # Number of records following header.
    N = dict(zones=14, boxes=80)[content]
    i = None
    for i,r in enumerate(f):
        data = struct.unpack(descriptor, r)
        suffix = data[-1]
        if format == 'v2':
            if 'zones' == content:
                title = id11fromzone(suffix)
            else:
                title = id11frombox(suffix, bos=bos)
        else:
            title = suffix
            output.write(title + '\n')
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
    if i is None:
        raise Error('No records found in file %r.' % inp.name)
    if i < N-1:
        raise Error('Unexpected end of file.')
    if i >= N:
        raise Error('Too many records.')

def id11fromzone(s):
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

def id11frombox(s, bos):
    """Convert a record suffix (5 words) as it appears in the BX file
    into an 11 character V2 Mean style station identifier.

    The returned strings will be of the form BOX@+NN+EEE where "+NN" and
    "+EEE" gives the lat/lon of the box centre.
    """

    d = struct.unpack("%s5i" % bos, s)
    box = d[1:5]
    centre = eqarea.centre(box)
    return "BOX@%+03.0f%+04.0f" % centre

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

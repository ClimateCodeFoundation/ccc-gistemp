#!/usr/bin/env python
# $URL$
# $Rev$
#
# subboxtopng.py
#
# David Jones, Clear Climate Code, 2010-08-26
#
# REQUIREMENTS
#
# pypng, from http://code.google.com/p/pypng/

"""
Convert subbox files to PNG image files.

subboxpng.py [--date YYYY-MM] [mask-or-subbox]

Converts either a mask file (text format, see work/step5mask for
example) or a subbox file (binary format,
result/SBBX1880.Ts.GHCN.CL.PA.1200 for example) into a PNG file.

Specifying a month with the --date option means that a subbox file is
converted, the PNG represents the specified date.  Otherwise a mask is
converted.

When converting a step5mask file the output is white for 0.000 (no land)
and black for 1.000 (use land).
"""

# http://docs.python.org/release/2.4.4/lib/module-itertools.html
import itertools

# http://code.google.com/p/pypng/
import png

# Clear Climate Code
import extend_path
from code import eqarea
from code.giss_data import MISSING

def topng(inp, date=None):
    """Convert *inp* into a PNG file.  Input file can be a step5mask
    file (produces greyscale PNG), or if *date* is supplied it can be a
    subbox file."""

    # :todo: move into proper module.
    from landmask import centrein

    resolution = 0.25
    width = 360/resolution
    height = 180/resolution
    assert int(width) == width
    assert int(height) == height
    width = int(width)
    height = int(height)

    cells = eqarea.grid8k()
    if not date:
        # Mask file in text format.
        values = maskbybox(inp, cells)
        colour = greyscale
        isgrey = True
    else:
        values = extractdate(inp, cells, date)
        colour = colourscale
        isgrey = False

    if isgrey:
        black = 0
    else:
        black = (0,0,0)
    # an array of rows:
    a = [[black]*width for _ in range(height)]

    for v,box in values:
        v = colour(v)
        for x,y in centrein(box, resolution):
            a[y][x] = v

    # For colour images each row of *a* is of the form:
    # [(R,G,B), (R,G,B), ...] we want to flatten it to:
    # [R,G,B,R,G,B,...]
    a = [list(itertools.chain(*row)) for row in a]

    try:
        outpath = inp.name + '.png'
    except:
        outpath = 'out.png'
    w = png.Writer(width=width, height=height,
      greyscale=isgrey, alpha=False,
      bitdepth=8)
    w.write(open(outpath, 'wb'), a)

def maskbybox(inp, grid):
    """Read a step5mask file box by box.  Yield (value, box) pair.
    """
    for row,box in itertools.izip(inp, grid):
        lat = float(row[:5])
        lon = float(row[5:11])
        s,n,w,e = box
        # If either of these fail, the input mask is in wrong sequence.
        assert s < lat < n
        assert w < lon < e
        v = float(row[16:21])
        yield v, box

def greyscale(v):
    """Convert value *v* in range 0 to 1 to a greyscale.  0 is white.
    """

    return 255-int(round(v*255))

def extractdate(inp, cells, date):
    """*date* should be a string in ISO 8601 format: 'YYYY-MM'.  From
    the binary subbox file *inp* extract the values corresponding to the
    date box by box.
    """
    import gio

    year,month = map(int, date.split('-'))

    records = iter(gio.SubboxReader(inp))
    meta = records.next()
    base_year = meta.yrbeg
    # Index of required month in the record series.
    i = (year - base_year)*12 + month - 1
    for record,box in itertools.izip(records, cells):
        assert record.first_year == base_year
        if i >= len(record.series):
            yield MISSING, box
        else:
            yield record.series[i], box

def colourscale(v):
    """Convert value *v* to a colour scale."""

    scale = [(-4, (0,0,255)), (0, (255,255,255)), (4, (255,0,0))]

    if v == MISSING:
        return (128,128,128)

    for i,knot in enumerate(scale):
        if knot[0] >= v:
            break
    else:
        i = len(scale)
    if i >= len(scale):
        return scale[-1][1]
    if i == 0:
        return scale[0][1]

    # *l* and *r* are the left and right knot values and colours.
    l = scale[i-1]
    r = scale[i]
    # Compute *t* the fractional distance between l[0] and r[0] (at
    # which *v* lies).
    t = (v-l[0]) / (r[0]-l[0])
    # Compute the interpolated colour.
    c = map(lambda a,b: lerp(t, a, b), l[1], r[1])
    c = [int(round(x)) for x in c]
    return c

def lerp(t, a, b):
    """Interpolate between *a* and *b*.  Result is *a* when *t* is 0;
    result is *b* when *t* is 1.
    """

    return (1-t)*a + t*b

def main(argv=None):
    import sys
    import getopt

    if argv is None:
        argv = sys.argv

    k = {}
    opt, arg = getopt.getopt(argv[1:], '', ['date='])
    for o,v in opt:
        if o == '--date':
            k['date'] = v
    if arg:
        for p in arg:
            topng(open(p), **k)
    else:
        topng(sys.stdin)

if __name__ == '__main__':
    main()

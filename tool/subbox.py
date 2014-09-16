#!/usr/bin/env python
#
# subbox.py
#
# David Jones, Clear Climate Code, 2010-08-26
# David Jones, Climate Code Foundation, 2014-09-15
#
# REQUIREMENTS
#
# pypng, from https://github.com/drj11/pypng

"""
Map a subbox file. Converts to rectangular PNG image; or,
when --polar is used, to a polar SVG image.

subbox.py [--polar] [--date YYYY-MM] [mask-or-subbox]

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

import math

# Clear Climate Code
import extend_path
from code import eqarea
from code.giss_data import MISSING
import gio

def to_polar_svg(inp, date=None):
    project = polar_project

    values = cells(inp, date)

    if not date:
        colour = greyscale
    else:
        colour = colourscale

    print """<svg
  xmlns="http://www.w3.org/2000/svg"
  xmlns:xlink="http://www.w3.org/1999/xlink"
  version="1.1">
"""
    print """<g transform="translate(270 270)">"""
    for v, rect in values:
        (s, n, w, e) = rect
        # :todo:(drj) For now, Northern Hemisphere only.
        if s < 0:
            continue
        ps = [(n, w), (n, e), (s, e), (s, w)]
        qs = [project(p) for p in ps]
        cell_svg(qs, colour(v))
    print "</g></svg>"

def polar_project(p):
    lat, lon = [math.radians(c) for c in p]
    # z = math.sin(lat)
    r = math.cos(lat)
    x = math.cos(lon) * r
    y = math.sin(lon) * r
    return x, y

def cell_svg(qs, fill_arg):
    """
    The argument qs is a list of 4 corners of a latitudinally
    oriented rectangular cell. Each corner is an (x,y)
    coordinate on a unit disc.

    Output a fragment of SVG.
    """

    import math

    s = 250

    # Radius of Northern edge.
    rn = math.hypot(*qs[0])
    # Radius of Southern edge.
    rs = math.hypot(*qs[-1])

    rn *= s
    rs *= s
    qs = [(x*s, -y*s) for x,y in qs]

    d = "M {:.1f} {:.1f} A {:.1f} {:.1f} {} {} {} {:.1f} {:.1f} L {:.1f} {:.1f} A {:.1f} {:.1f} {} {} {} {:.1f} {:.1f} z".format(
      qs[0][0], qs[0][1],
      rn, rn, 0, 0, 1, qs[1][0], qs[1][1],
      qs[2][0], qs[2][1],
      rs, rs, 0, 0, 1, qs[3][0], qs[3][1])

    fill = "#{:02x}{:02x}{:02x}".format(*fill_arg)

    print """<path fill="{}" stroke="none" d='{}' />""".format(fill, d)


def to_rect_png(inp, date=None):
    """
    Convert *inp* into a PNG file.  Input file can be a step5mask
    file (produces greyscale PNG), or if *date* is supplied it can be a
    subbox file.
    """

    # http://code.google.com/p/pypng/
    import png

    values = cells(inp, date)

    # :todo: move into proper module.
    from landmask import centrein

    resolution = 0.25
    width = 360/resolution
    height = 180/resolution
    assert int(width) == width
    assert int(height) == height
    width = int(width)
    height = int(height)

    if not date:
        colour = greyscale
    else:
        colour = colourscale

    if colour == greyscale:
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
    if colour != greyscale:
        a = [list(itertools.chain(*row)) for row in a]

    try:
        outpath = inp.name + '.png'
    except:
        outpath = 'out.png'
    w = png.Writer(width=width, height=height,
      greyscale=(colour == greyscale), alpha=False,
      bitdepth=8)
    w.write(open(outpath, 'wb'), a)

def cells(inp, date=None):
    """
    Yield a series of (value, rect) pairs.
    """
    subboxes = eqarea.grid8k()
    if not date:
        # Mask file in text format.
        values = gio.maskboxes(inp, subboxes)
    else:
        values = extractdate(inp, subboxes, date)

    return values

def greyscale(v):
    """
    Convert value *v* in range 0 to 1 to a greyscale.  0 is white.
    """

    return 255-int(round(v*255))

def extractdate(inp, cells, date):
    """
    *date* should be a string in ISO 8601 format: 'YYYY-MM'.  From
    the binary subbox file *inp* extract the values corresponding to the
    date box by box.
    """

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
    """
    Convert value *v* to a colour scale.
    """

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
    """
    Interpolate between *a* and *b*.  Result is *a* when *t* is 0;
    result is *b* when *t* is 1.
    """

    return (1-t)*a + t*b

def to_ghcnm(inp):
    """
    Convert a file inp from subbox to GHCN-M (v3) format.
    """

    import sys

    out = sys.stdout

    w = gio.GHCNV3Writer(file=out)

    subbox = iter(gio.SubboxReader(inp))
    # First record is metadata, which we ignore.
    subbox.next()
    for record in subbox:
        lat,lon = eqarea.centre(record.box)
        record.uid = '%+05.1f%+06.1f' % (lat,lon)
        w.write(record)
    w.close()

def main(argv=None):
    import sys
    import getopt

    if argv is None:
        argv = sys.argv

    command = to_rect_png

    k = {}
    opt, arg = getopt.getopt(argv[1:], '', ['ghcnm', 'date=', 'polar'])
    for o,v in opt:
        if o == '--date':
            k['date'] = v
        if o == '--polar':
            command = to_polar_svg
        if o == '--ghcnm':
            command = to_ghcnm


    if arg:
        arg = [open(p) for p in arg]
    else:
        arg = [sys.stdin]

    for a in arg:
        command(a, **k)

if __name__ == '__main__':
    main()

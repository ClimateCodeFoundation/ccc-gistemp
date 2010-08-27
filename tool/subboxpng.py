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

Conceptually intended to convert a variety of subbox files, but for now
converts step5mask files.

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

def topng(inp):
    """Convert a step5mask into a PNG file."""

    # :todo: move into proper module.
    from landmask import centrein

    resolution = 0.25
    width = 360/resolution
    height = 180/resolution
    assert int(width) == width
    assert int(height) == height
    width = int(width)
    height = int(height)
    # an array of rows:
    a = [[0]*width for _ in range(height)]

    for row,box in itertools.izip(inp, eqarea.grid8k()):
        lat = float(row[:5])
        lon = float(row[5:11])
        s,n,w,e = box
        # If either of these fail, the input mask is in wrong sequence.
        assert s < lat < n
        assert w < lon < e
        v = float(row[16:21])
        v = 255-int(round(v*255))
        for x,y in centrein(box, resolution):
            a[y][x] = v

    try:
        outpath = inp.name + '.png'
    except:
        outpath = 'out.png'
    w = png.Writer(width=width, height=height,
      greyscale=True, alpha=False,
      bitdepth=8)
    w.write(open(outpath, 'wb'), a)

def main(argv=None):
    import sys
    # http://docs.python.org/release/2.4.4/lib/module-fileinput.html
    import fileinput

    if argv is None:
        argv = sys.argv

    arg = argv[1:]
    if arg:
        for p in arg:
            topng(open(p))
    else:
        topng(sys.stdin)

if __name__ == '__main__':
    main()

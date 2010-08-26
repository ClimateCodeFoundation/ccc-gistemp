#!/usr/bin/env python
# $URL$
# $Rev$
#
# maskpng.py
#
# David Jones, Clear Climate Code, 2010-08-26

"""
maskpng.py land_percent.asc

Convert the land percentage grid to PNG image.

The input file should be an ASCII text land percentage file.
First row is northernmost with cells running from -180
to +180.  Resolution of file is computed from file format
(360 values per row is 1 degree, 720 is 0.5 degree, 1440 is 0.25 degree)

REFERENCES

http://islscp2.sesda.com/ISLSCP2_1/data/ancillary/land_water_masks_xdeg/0_land_water_masks_readme.txt
http://islscp2.sesda.com/ISLSCP2_1/html_pages/groups/ancillary/land_water_masks_xdeg.html

"""

import math

# http://code.google.com/p/pypng/
import png

def topng(inp, out):
    """Given a land percent file as input, produce a PNG image
    as output."""

    land = grid(inp)
    w = png.Writer(width=land.w, height=land.h,
      greyscale=True, alpha=False,
      bitdepth=8)
    w.write(out, ([int(round(x/100.0*255)) for x in row] for row in land.a))

# Copied from landmask.py
def grid(inp):
    """Convert ASCII file to regular grid."""

    class Struct:
        pass
    gridded = Struct()
    gridded.a = [map(int, row.split()) for row in inp]
    gridded.w = len(gridded.a[0])
    gridded.h = len(gridded.a)
    gridded.resolution = 360.0/gridded.w
    assert gridded.resolution in (1.0, 0.5, 0.25)
    return gridded

def main(argv=None):
    import sys
    if argv is None:
        argv = sys.argv

    arg = argv[1:]
    if not arg:
        print __doc__
        return 2

    topng(open(arg[0]), sys.stdout)

if __name__ == '__main__':
    main()

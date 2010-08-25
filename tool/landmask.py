#!/usr/bin/env python
# $URL$
# $Rev$
#
# landmask.py
#
# David Jones, Clear Climate Code, 2010-08-25

"""
landmask.py land_percent.asc

landmask computes an 8000 cell land mask suitable for using as the input
to ccc-gistemp Step 5.  The input file should be an ASCII text land
percentage file.  First row is northernmost with cells running from -180
to +180.  Resolution of file is computed from file format
(360 values per row is 1 degree, 720 is 0.5 degree, 1440 is 0.25 degree)

REFERENCES

http://islscp2.sesda.com/ISLSCP2_1/data/ancillary/land_water_masks_xdeg/0_land_water_masks_readme.txt
http://islscp2.sesda.com/ISLSCP2_1/html_pages/groups/ancillary/land_water_masks_xdeg.html

"""

import math

# Clear Climate Code
import extend_path
from code import eqarea
from code import giss_data

def maskit(inp, out):
    """Given a land percent file as input, produce a GISTEMP cell mask
    as output."""

    land = grid(inp)
    resolution = land.resolution
    for subbox in eqarea.grid8k():
        values = [land.a[y][x] for x,y in centrein(subbox, resolution)]
        # For GISTEMP we mask as land if there is _any_ land in the
        # cell.
        if sum(values) > 0:
            mask = 1
        else:
            mask = 0
        out.write("%sMASK%.3f\n" % (giss_data.boxuid(subbox), mask))

def centrein(box, resolution):
    """Grid a sphere with cells spaced every *resolution* degrees in
    latitude and longitude, then return the integer coordinates of those
    cells that have a centre that lies in *box*.
    
    Coordinates are returned as an (x,y) pair where (0,0) is the
    northernmost, westernmost (being -180) cell, and *x* corresponds
    to latitude.
    """

    s,n,w,e = box
    for y in range(math.floor((s + 90.0)/resolution + 0.5),
                   math.floor((n + 90.0)/resolution + 0.5)):
        for x in range(math.floor((w + 180.0)/resolution + 0.5),
                       math.floor((e + 180.0)/resolution + 0.5)):
            # Flip y
            yield x,int(90/resolution - y - 1)

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

    maskit(open(arg[0]), sys.stdout)

if __name__ == '__main__':
    main()

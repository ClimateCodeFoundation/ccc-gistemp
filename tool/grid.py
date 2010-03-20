#!/usr/bin/env python
# $URL$
# $Rev$
#
# grid.py

"""
grid YYYY-MM [v2-file]

Display gridded anomalies as SVG file.
"""

def map(when, inp, out):
    """Take a ccc-gistemp subbox file in V2 mean format as *inp* and
    produce an SVG file on *out*."""

    import math
    import re

    out.write("""<svg 
      xmlns="http://www.w3.org/2000/svg"
      xmlns:xlink="http://www.w3.org/1999/xlink"
      version="1.1">\n""")

    m = re.match(r'(\d{4})-(\d{2})', when)
    year, month = m.groups()
    # Convert to 0-based month.
    month = int(month) - 1
    assert 0 <= month < 12

    out.write("""<g transform='translate(180,90)'>\n""")
    for (lat,lon),v in filter_month(inp, year, month):
        x,y = topixel(lat, lon)
        y = -y
        if v > 0:
            fill = 'red'
        else:
            fill = 'blue'
        radius_scale = 0.25
        r = math.sqrt(abs(v)) * radius_scale
        out.write("""<circle cx='%.1f' cy='%.1f' r='%.1f' fill='%s' />\n""" %
          (x, y, r, fill))
    out.write('</g>\n')
    out.write("""</svg>\n""")

def topixel(lat, lon):
    u"""Return x,y coordinate.  Plate Carr\xe9e projection."""

    return lon,lat


def filter_month(inp, year, month):
    """Yield each value from the v2.mean file for the year and month in
    question.  Invalid values (marked with -9999 in the file) are
    ignored.  Each value is yielded as ((lat,lon),v) where lat,lon are
    the latitude and longitude in degrees (positive is North and East).
    """

    for line in inp:
        if line[12:16] == year:
            v = int(line[16+5*month:21+5*month])
            if v != -9999:
                lat = float(line[:5])
                lon = float(line[5:11])
                yield ((lat,lon), v)

def main(argv=None):
    import sys
    if argv is None:
        argv = sys.argv

    when = argv[1]
    if len(argv) > 2:
        f = open(argv[2])
    else:
        f = sys.stdin
    map(when, f, sys.stdout)

if __name__ == '__main__':
    main()

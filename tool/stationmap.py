#!/usr/bin/env python
# $URL$
# $Rev$
# stationmap.py
# David Jones, 2011-08-13, Climate Code Foundation.

"""
stationmap.py station-list

Plot the stations on a map (SVG output).

*station-list* names a file of station identifiers.  The file is scanned
and any line starting with 11 alphanumeric characters (no spaces) is
used.  The location metadata for the stations is taken from the file
input/v2.inv.
"""

import re
import sys

# ccc-gistemp
import gio

def map(inp, out=sys.stdout):
    """Plot map.  *inp* is the name of the file of station IDs."""

    out.write("""<svg
      xmlns="http://www.w3.org/2000/svg"
      xmlns:xlink="http://www.w3.org/1999/xlink"
      version="1.1">
""")
    out.write("""<defs>
  <style type="text/css">
    g.stations {
      stroke: green; stroke-width: 1.4; stroke-linecap: round; fill: none }
    text { fill: black; font-family: Verdana }
  </style></defs>
""")

    out.write("<g class='stations'>\n")
    out.write("<g transform='translate(0, 360) scale(1,-1)'>\n")
    for station in stations(inp):
        lat,lon = station.lat,station.lon
        lon += 180.0
        lat += 90.0
        lat,lon = [x*2.0 for x in (lat,lon)]
        out.write("<path d='M%.2f %.2fl0 0' />\n" % (lon,lat))
    out.write("</g>\n"*2)
    out.write("</svg>\n")

def stations(filenames):
    """Yield the metadata for each station listed in *filenames*."""

    filename = filenames[0]

    ids = set(l[:11] for l in open(filename) if re.match(r'\w{11}', l))
    meta = gio.station_metadata(path='input/v2.inv')
    lost = False
    for id in ids:
        if id in meta:
            yield meta[id]
        else:
            lost = True
    if lost:
        raise Exception("Couldn't find metadata for some stations.")


def main(argv=None):
    import sys
    if argv is None:
        argv = sys.argv

    arg = argv[1:]

    map(arg)

if __name__ == '__main__':
    main()

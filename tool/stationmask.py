#!/usr/bin/env python
# $URL$
# $Rev$
#
# stationmask.py
#
# David Jones, Clear Climate Code, 2010-08-27

"""
stationmask.py cells [stations]

Extract a list of stations that are in a selected set of cells.  The subbox
cells are taken from the *cells* file which should be in the same
format as step5mask (value 0.000 excludes that cell, any other value
includes it).  For rectangular regions, you can produce a suitable mask
file with tool/rectmask.py.

The locations of the stations are taken from the input/v2.inv file.  The
list of stations is taken from the *stations* file, using input/v2.inv
when this is not specified.
"""

# Clear Climate Code
import extend_path
from code import eqarea

def stationmask(inp, stations, out, inv, dribble):
    """Take list of cells from *inp*, list of stations from *stations*,
    and output on *out* a list of station identifiers (11-digit) (that
    are within the cells specified in *inp); locations are taken from
    the v2.inv file *inv*.
    """

    centres = [row[:11] for row in inp if float(row[16:21])]
    print >> dribble, len(centres), "centres"
    centres = [(float(c[:5]),float(c[5:11])) for c in centres]
    boxes = eqarea.grid8k()
    selectboxes = [box for box in boxes if
      [centre for centre in centres if boxcontains(box, centre)]]
    print >> dribble, len(selectboxes), "boxes"
    # Set of stations.
    stations = set(l[:11] for l in stations)
    meta = list(inv)
    selected = []
    # Iterate over all entries in v2.inv, discarding those not in
    # *stations*.
    for line in meta:
        uid = line[:11]
        if uid not in stations:
            continue
        lat = float(line[43:49])
        lon = float(line[50:57])
        for box in selectboxes:
            if eqarea.boxcontains(box, (lat,lon)):
                selected.append(uid)
                dribble.write('\r%d stations' % len(selected))
                break
    dribble.write('\n')
    selected = set(selected)
    for id11 in sorted(selected):
        out.write(id11 + '\n')

def usage():
    print __doc__
    return 2

def main(argv=None):
    import sys
    if argv is None:
        argv = sys.argv
    arg = argv[1:]
    if len(arg) < 1:
        return usage()
    if len(arg) > 1:
        stations = arg[1]
    else:
        stations = 'input/v2.inv'

    stationmask(open(arg[0]), stations=open(stations),
      out=sys.stdout, inv=open('input/v2.inv'),
      dribble=sys.stderr)

if __name__ == '__main__':
    main()

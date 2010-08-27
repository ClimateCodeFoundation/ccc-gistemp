#!/usr/bin/env python
# $URL$
# $Rev$
#
# stationmask.py
#
# David Jones, Clear Climate Code, 2010-08-27

"""
stationmask.py cells

From the input/v2.inv file extract a list of stations that are in a
selected set of cells.  The subbox cells are taken from the cells file
which should be in the same format as step5mask (value 0.000 excludes
that cell, any other value includes it).
"""

# Clear Climate Code
import extend_path
from code import eqarea

def stationmask(inp, out, inv, dribble):
    """Take list of cells from *inp*, and output on *out* a list of
    station identifiers (11-digit), taken from the v2.inv file *inv*.
    """

    centres = [row[:11] for row in inp if float(row[16:21])]
    print >> dribble, len(centres), "centres"
    centres = [(float(c[:5]),float(c[5:11])) for c in centres]
    boxes = eqarea.grid8k()
    selectboxes = [box for box in boxes if
      [centre for centre in centres if boxcontains(box, centre)]]
    print >> dribble, len(selectboxes), "boxes"
    stations = list(inv)
    selected = []
    for line in stations:
        lat = float(line[43:49])
        lon = float(line[50:57])
        for box in selectboxes:
            if boxcontains(box, (lat,lon)):
                selected.append(line[:11])
                dribble.write('\r%d stations' % len(selected))
                break
    dribble.write('\n')
    selected = set(selected)
    for id11 in sorted(selected):
        out.write(id11 + '\n')

def boxcontains(box, p):
    """True iff *box* (4-tuple of (s,n,w,e) ) contains point *p* (pair
    of (lat,lon)."""

    s,n,w,e = box
    return s <= p[0] < n and w <= p[1] < e

def main(argv=None):
    import sys
    if argv is None:
        argv = sys.argv
    arg = argv[1:]

    stationmask(open(arg[0]), sys.stdout, open('input/v2.inv'),
      sys.stderr)

if __name__ == '__main__':
    main()

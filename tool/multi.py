#!/usr/bin/env python
# $URL$
# $Rev$

"""
Multi-tool.  A collection of miscellaneous commands and procedures.
Mostly intended for debugging and introspection.

Code is collected here to avoid cluttering up tool with many one-off
scripts.

Run "python tool/multi.py commands" to see command list.
"""

import extend_path

# http://docs.python.org/release/2.6.6/library/json.html
import json
# http://docs.python.org/release/2.4.4/lib/module-itertools.html
import itertools
# http://docs.python.org/release/2.4.4/lib/module-os.html
import os
import re
import sys
import urllib

# Clear Climate Code
import extend_path
import gio
from code.giss_data import valid, MISSING

class Fatal(Exception):
    """An error occurred."""

def tabletov2(argv):
    """[command] Convert GISTEMP table format (GLB.* and so on) to GHCN
    v2 format.
    """

    dotabletov2(arg=argv[1:], out=sys.stdout)

def dotabletov2(arg, out, labels=None):
    """Open the files named in the list *arg* and write out GHCN v2
    format to *out*.  If *labels* is specified, it should be a list of
    12 character strings to be used for the station identifiers in the
    GHCN v2 output file; one element per input file.
    """

    v2 = gio.GHCNV2Writer(file=out, scale=1)
    if not labels:
        labels = [name[:12] for name in arg]
    for name,uid in zip(arg,labels):
        f = urllib.urlopen(name)
        for year,monthlies in tablemonthlies(f):
            v2.writeyear(uid, year, monthlies)
    v2.close()

def tablemonthlies(f):
    """Convert GISTEMP table file (GLB.*) to sequence of (year, seq)
    pairs.  Each *seq* is a list of 12 monthly values.
    """

    def convert(s):
        try:
            return int(s)
        except:
            return MISSING

    for row in f:
        if not re.match(r'\d{4}', row):
            continue
        year = int(row[:4])
        monthlies = [convert(row[i:i+5]) for i in range(5, 65, 5)]
        yield year,monthlies

def timecount(argv):
    """[command] [--mask maskfile] [--dir .] Count of stations
    used over time."""

    # http://docs.python.org/release/2.4.4/lib/module-getopt.html
    import getopt
    opts,arg = getopt.getopt(argv[1:], '', ['dir=', 'mask='])
    option = dict((o[2:],v) for o,v in opts)

    counts(out=sys.stdout, **option)

def counts(out, dir='.', mask=None):
    """Print to *out* a count for each month of the number of stations
    used.
    
    *dir* specifies the name of a directory for the input files
    (that are intermediate products of ccc-gistemp);
    from this directory the files used are: 'log/step3.log' (to
    determine what stations and months are used), and 'work/step2.v2'
    (to determine which specific months of those stations have data
    present).

    *mask* specifies the name of a mask file to be used; only the
    cells present in the mask file will be examined for station usage.
    """

    # Dictionary that maps from station (12 char identifier) to
    # monthset, where *monthset* is a set of the numbers from 0
    # (January) to 11 (December).
    stationmonths = dict()

    log = os.path.join(dir, 'log', 'step3.log')
    step2 = os.path.join(dir, 'work', 'step2.v2')

    for row in stations_logged(mask, log=log):
        stations = json.loads(row[2])
        for station,weight,months in stations:
            stationmonths[station] = monthset(months) | stationmonths.get(
              station, set())

    monthcount = {}
    for station in gio.GHCNV2Reader(path=step2):
        if station.uid not in stationmonths:
            continue
        for m in stationvalidmonths(station):
            monthcount[m] = monthcount.get(m, 0) + 1
    v2monthcount(out, monthcount)

def v2monthcount(out, counts):
    """Write *counts* to *out* in GHCN V2 format.  *counts* is a dict
    that maps from *midx* to count, whre *midx* is a 0 based index for
    the months (year 1 January is 12).
    """

    ymin = min(counts)//12
    ymax = max(counts)//12
    id12 = 'STATIONCOUNT'
    for y in range(ymin, ymax+1):
        vs = [counts.get(y*12+m, 0) for m in range(12)]
        s = ("%5d"*12) % tuple(vs)
        out.write("%s%04d%s\n" % (id12, y, s))

def monthset(s):
    """Convert string of form "001111110000" to a set of months.  For
    each character in the string that is '1', the corresponding index is
    added to the result set.
    """

    assert len(s) == 12
    assert set(s) <= set('01')

    return set(i for i,x in enumerate(s) if x == '1')

def stationvalidmonths(record):
    """Return the set of months for which the record has valid data.
    Each month is encoded as a number with january of year 1 being 12."""

    first = record.first_month - 1
    return set(first+i for i,v in enumerate(record.series) if valid(v))

def wherestations(argv):
    """[command] [--inv v2.inv] stations [...]  Outputs the cells
    (as a mask file) that are occupied by at least one station list in
    any of the *stations* files.
    """

    # http://docs.python.org/release/2.4.4/lib/module-getopt.html
    import getopt
    opts,arg = getopt.getopt(argv[1:], '', ['inv='])
    option = dict((o[2:],v) for o,v in opts)
    where_stations(arg, out=sys.stdout, **option)

here = os.path.abspath(__file__)
parent = os.path.dirname(os.path.dirname(here))
input = os.path.join(parent, 'input')

def where_stations(stations=[], out=None,
  inv=os.path.join(input, 'v2.inv')):
    """Implementation of whatstations command."""

    from code import eqarea
    from code import giss_data

    stations = set(l[:11] for l in
      itertools.chain(*[open(f) for f in stations]) if ':' not in l)
    meta = gio.station_metadata(path=inv)
    assert meta
    # Restrict meta to the set *stations*
    meta = dict((uid,m) for uid,m in meta.iteritems() if uid in stations)

    count = eqarea.GridCounter()
    for s in meta.itervalues():
        count(s.lat, s.lon)
    for c,cell in count.boxes():
        m = c > 0
        out.write("%sMASK%.3f\n" % (giss_data.boxuid(cell), m))

def whatstations(argv):
    """[command] [--mask maskfile] [--log step3.log] Determines
    what (land) stations are used for a particular area.

    whatstations --mask maskfile

    Works by examining the Step 3 log file, which can be specified with
    --log option.
    """

    # http://docs.python.org/release/2.4.4/lib/module-getopt.html
    import getopt
    opts,arg = getopt.getopt(argv[1:], '', ['log=', 'mask='])
    option = dict((o[2:],v) for o,v in opts)
    stations(out=sys.stdout, **option)

def stations_logged(mask=None, log='log/step3.log'):
    """An iterator that yields each of the log records from step3 that
    identify the stations used (that is, the second item in the row is
    "stations").  *mask* (optionally) specifies the name of a cell mask
    file; when used only the cells present in the mask are examined.
    *log* specifies the name of the log file to examine.

    Each log record is yielded as a triple (cellid, 'stations',
    JSONstring).  *JSONstring* is a string suitable for passing to
    json.loads().
    """

    log = open(log)
    if mask:
        cells = cellsofmask(open(mask))

    for row in log:
        row = row.split(' ', 2)
        if row[1] != 'stations':
            continue
        if not mask or row[0] in cells:
            yield row

def stations(out, log='log/step3.log', mask=None):
    """Print to *out* a list of the stations used."""

    station_weight = dict()
    station_months = dict()
    cellcount = 0
    for row in stations_logged(log=log, mask=mask):
        stations = json.loads(row[2])
        for item in stations:
            station,weight,months = item[:3]
            if weight:
                station_weight[station] = max(
                  station_weight.get(station, 0), weight)
                station_months[station] = station_months.get(
                  station, set()) | monthset(months)
        cellcount += 1
    
    def asstr(s):
        """Convert set of months to string."""

        return ''.join('01'[i in s] for i in range(12))

    if mask:
        ncells = len(cellsofmask(open(mask)))
    else:
        ncells = 8000
    print >>out, "Cells: %d/%d" % (cellcount, ncells)
    print >>out, "Stations: %d" % len(station_weight)
    for station,weight in sorted(station_weight.iteritems()):
        print >>out, station, asstr(station_months[station]), weight

def cellsofmask(inp):
    """*inp* is an open mask file.  The returned value is a set of all
    the cells (identified by their 12-character ID) that have a non-zero
    value."""

    return set(row[:12] for row in inp if float(row[16:21]))

def maskmap(arg):
    """[command] mask  Draws an SVG map of the cells in mask."""

    from code import eqarea

    maskfile = arg[1]
    out = sys.stdout

    def colour(x):
        """Choose colour for cell.  Returned as an (r,g,b,a) 4-tuple
        with values between 0.0 and 1.0.
        """

        assert v in (0,1)
        if v:
            return (1.0, 0.0, 0.0, 0.4)
        else:
            return (0.0, 0.0, 0.0, 0.0)
    def SVGfill(t):
        """Render the 4-tuple t as an attribute fragment that specifies
        fill and fill-opacity.  For example,
        "fill='rgb(192,192,255)' fill-opacity='1.0'".
        """

        return ("fill='rgb(%.0f,%.0f,%.0f)' " % tuple(255.0*x for x in t[:3]) +
          "fill-opacity='%.3f'" % t[3])

    def transform(p):
        u"""Transform p=(lat,lon) into (x,y) used for map system.  Which
        in this case is Plate Carr\xe9e with x from 0 to 720, and y from
        0 to 360."""

        lat,lon = p
        y = (90+lat)*2
        x = (180+lon)*2
        return x,y

    out.write("""<svg 
      xmlns="http://www.w3.org/2000/svg"
      xmlns:xlink="http://www.w3.org/1999/xlink"
      version="1.1">\n""")

    out.write("<g class='mask' transform='translate(0,360) scale(1,-1)'>\n")
    cells = eqarea.grid8k()
    for v,box in gio.maskboxes(open(maskfile), cells):
        s,n,w,e = box
        corners = [(n,w), (n,e), (s,e), (s,w)]
        corners = map(transform, corners)
        out.write("<path d='M %.2f %.2f L" % corners[0] +
          " %.2f"*6 % tuple(itertools.chain(*corners[1:])) +
          " z' %s />\n" % SVGfill(colour(v)))
    out.write("</g>\n")
    out.write("</svg>\n")


def cmpcontributors(arg):
    """[command] Compares the lists of stations for each cell (if using
    step3.log files); or the lists of cells for each region (if using
    step5.log files).  arg[1] and arg[2] should be the names of two
    different log files.
    """

    fs = map(open, arg[1:3])
    stations = map(celldict, fs)
    if set(stations[0]) != set(stations[1]):
        print "Sets of cells differ"
    common = set(stations[0]) & set(stations[1])
    print "%d common cells" % len(common)
    for cell in common:
        # Get the two collections of stations...
        # as lists...
        stationsa = stations[0][cell]
        stationsb = stations[1][cell]
        # as dicts...
        dicta = dict(stationsa)
        dictb = dict(stationsb)
        # and sets.
        seta = set(dicta)
        setb = set(dictb)
        if stationsa != stationsb:
            ina = seta - setb
            inb = setb - seta
            reportinonelist(cell, fs[0].name, ina, dicta)
            reportinonelist(cell, fs[1].name, inb, dicta)
            commonstations = seta & setb
            for station in commonstations:
                if dicta[station] != dictb[station]:
                    print "Cell %s station %s has different weights (%f and %f)." %(
                      cell, station, dicta[station], dictb[station])
            # Split into stations with 0 weight and those with non-zero
            # weight.
            zeroa = [s for s,w in stationsa if not w and s in commonstations]
            zerob = [s for s,w in stationsb if not w and s in commonstations]
            posa = [s for s,w in stationsa if w and s in commonstations]
            posb = [s for s,w in stationsb if w and s in commonstations]
            for i,(a,b) in enumerate(zip(posa, posb)):
                if a != b:
                    print "Cell %s order differs: %s and %s at index %d" % (
                       cell, a, b, i)
                break
            # Not normally reported.
            if False and zeroa != zerob:
                print "Cell %s (the unused stations are in different order; result not affected)" % (
                  cell)


def celldict(f):
    """From the (open) step3.log file *f* return a dict that maps
    from box identifier (12 characters) to a list of (station,weight)
    pairs.
    """

    result = {}
    for row in f:
        row = row.split(' ', 2)
        # The 2nd item is 'stations' for step3.log and 'cells' for
        # step5.log.
        if row[1] not in ['stations', 'cells']:
            continue
        stations = json.loads(row[2])
        pairs = [(t[0],t[1]) for t in stations]
        result[row[0]] = pairs
    return result

def reportinonelist(cell, name, culprit, weight):
    """For the cell and logfile identified by the strings *cell*
    and *name*, print a report of the stations contributing to
    that cell that appear only in that file (and not in the
    corresponding cell in the other file).  *culprit* is a
    sequence of the station identifiers, *weight* is a dict that
    maps from station identifiers to weight.
    """

    if not culprit:
        return
    # Split culprits into those with zero weight, and those with
    # non-zero weight
    zerow = [c for c in culprit if weight[c] == 0]
    somew = [c for c in culprit if weight[c] != 0]
    if zerow:
        print "Cell %s %s only (zero weights): %s " % (cell, name,
          ' '.join(sorted(zerow)))
    if somew:
        print "Cell %s %s only: %s " % (cell, name,
          ' '.join(sorted(somew)))


def blog20110106(arg):
    """Generates KML file for blog post.  One of the ccc-gistemp boxes,
    and the extended box used to do station culling.
    """

    from code import earth
    import math
    from code import eqarea

    p = (64.2,-90)
    d = 1200.0/earth.radius
    dd = math.degrees(d)

    circle = surfcircle(p, d)
    circle.append(circle[0])

    allboxes = list(eqarea.gridsub())
    # Let *box* be Boundaries of box in order: S, N, W, E.
    box,cells = allboxes[6]
    # Let *cells* be all the cells in that box.
    # Get the Northwest most cell.
    nwcell = max(cells, key=lambda x:x[1]-x[2])
    nwcentre = eqarea.centre(nwcell)
    s,n,w,e = box
    # List of coordinates for box 7.
    box7 = [(n,w), (n,e), (s,e), (s,w)]
    box7.append(box7[0])

    # "Extended" box.
    span = e-w
    ebox = [s-dd, n+dd, w-span/2, e+span/2]
    s,n,w,e = ebox
    # List of coordinate for extended box 7.
    ebox7 = [(n,w), (n,e), (s,e), (s,w)]
    ebox7.append(ebox7[0])

    circle2 = surfcircle(nwcentre, d)
    circle2.append(circle2[0])

    print """<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
<Document>
  <name>Station Culling</name>
""" + kmlpolystyle(color='aa0000ff', id='s0') + """
""" + kmlpolystyle(color='aaffffff', id='s1') + """
""" + kmlpolystyle(color='5500ff00', id='s2') + """
""" + kmlplace("1200 km circle", circle, style="#s0") + """
""" + kmlplace("circle at %r" % (nwcentre,), circle2, style="#s1") + """
""" + kmlplace("Culling Box", ebox7, style="#s2") + """
""" + kmlplace("Box 7", box7, style="#s1") + """
</Document>
</kml>"""

def nature201102(arg):
    """Generate SVG for image as per Nature's request.
    (expects to find a completed result directory).
    """

    # Download GISTEMP global result and copy together with the global
    # result on local disk to a GHCN v2 format file, 'nature.v2'
    out = open('nature.v2', 'w')
    labels=['_ccc-gistemp', 'NASA_GISTEMP', 'cccNASA_diff']

    print "fetching data..."
    dotabletov2(['result/mixedGLB.Ts.ho2.GHCN.CL.PA.txt',
      'http://data.giss.nasa.gov/gistemp/tabledata/GLB.Ts+dSST.txt'], out,
      labels=labels[:2])
    out.close()

    # Compute difference series: ccc-gistemp - GISTEMP:
    import gio
    from code import giss_data

    reader = gio.GHCNV2Reader(path=out.name)
    stations = list(reader)
    assert 2 == len(stations)
    assert stations[0].first_year == stations[1].first_year
    difference = giss_data.Series()
    difference.set_series(stations[0].first_month,
      [p-q for p,q in zip(stations[0].series, stations[1].series)
        if p != MISSING and q != MISSING])
    difference.uid = labels[2]
    w = gio.GHCNV2Writer(file=open(out.name, 'a'))
    w.write(difference)
    w.close()

    print "generating SVG..."
    from tool import stationplot

    stationplot.main(
      ("stationplot -o nature.svg --axes=yyr --offset 0,-0.2,0 -c rscale=6000;yscale=300;ytick=0.2;legend=none;buginkscapepdf=1 -s 0.01 -y -d nature.v2 %s %s %s" % tuple(labels)).split())
    print "generating PDF..."
    # Expects to find "inkscape" on PATH.
    # On drj's Mac, add
    # "/Applications/Inkscape.app/Contents/Resources/bin" to the PATH
    # first.
    # See http://inkscape.modevia.com/inkscape-man.html for inkscape
    # command line use.
    os.system('inkscape --export-pdf=nature.pdf nature.svg')


def kmlpolystyle(color, id):
    """A simply Style including a PolyStyle with color."""

    assert 8 == len(color)

    return """<Style id='""" + id + """'>
    <PolyStyle>
      <color>""" + color + """</color>
    </PolyStyle>
  </Style>
  """

def kmlplace(name, coords, style=None):
    """Return a string suitable to be inserted into a kml document as a
    Placemark element.  *name* is a string; *coords* is a list of (lat,lon)
    coordinate pairs.  If *style* is specified, it should be the URL of
    the style element (usually "#something").
    """

    from xml.sax.saxutils import escape

    name = escape(name)
    coordstring = kmlcoord(coords)

    if style:
        styleurl = "<styleUrl>" + escape(style) + "</styleUrl>"
    else:
        styleurl = ""

    return """<Placemark>
  <name>""" + name + """</name>
""" + styleurl + """
    <Polygon>
      <outerBoundaryIs>
        <LinearRing>
          <coordinates>
""" + coordstring + """
          </coordinates>
        </LinearRing>
      </outerBoundaryIs>
    </Polygon>
  </Placemark>
"""


def kmlcoord(l):
    """Convert list of (lat,lon) pairs to a list suitable for the
    contents of a KML <coordinates> element.
    """

    return ' '.join(','.join(map(str, [p[1],p[0],0])) for p in l)

def surfcircle(point, arc):
    """Compute the circle that is *arc* radians from *point* (which is a
    (lat,lon) pair in degrees).  A list is returned which is a polygon
    that approximates the circle.  The returned list is a list of
    (lat,lon) pairs (in degrees); the elements of the list are distinct
    (in particular, the last element is not equal to the first).
    """

    import math

    # Proceed by:
    # - Choose some point *q*, where the angle between *point* and
    #     *q* is *arc* radians.  Convert *q* to 3D cartesian coords.
    # - Construct the 3x3 matrix that is a rotation about the vector
    #     through *point*, rotating through an angle of *s* degrees (s = 10,
    #     say).
    # - Repeatedly apply the matrix to *q* to get the orbit of *q* under
    #     the rotation.
    # - The list of points in the orbit is the required polygon; convert
    #   to (lat,lon) representation.

    s = 10

    p = [math.radians(x) for x in point]
    # If p is in Northern Hemisphere, pick q to be *arc* radians South;
    # and pick q to be *arc* radians North when p is in Southern
    # Hemisphere.
    if p[0] >= 0:
        q = (p[0]-arc, p[1])
    else:
        q = (p[0]+arc, p[1])

    pcart = cartesian3D(p)
    qcart = cartesian3D(q)

    s = 10
    th = math.radians(s)
    rot = m3x3rotVth(pcart, th)

    result = []
    for _ in range(0,360,s):
        result.append(qcart)
        qcart = m3x3multV(rot, qcart)

    return converttolatlon(result)

def converttolatlon(l):
    """Convert a list of (x,y,z) cartesian triples to list of (lat,lon)
    pairs (in degrees)."""

    return [convert1latlon(p) for p in l]

def convert1latlon(p):
    """Convert *p* from a cartesian triple (x,y,z) to (lat,lon) pair in
    degrees."""

    assert 0.99 < sum(v*v for v in p)**0.5 < 1.01

    import math

    x,y,z = p

    lon = math.atan2(y,x)
    lat = math.asin(z)

    return [math.degrees(v) for v in (lat,lon)]

def cartesian3D(p):
    """Convert *p*, which is a (lat,lon) pair in radians, to 3D
    cartesian coordinates."""

    import math

    lat,lon = p
    x = math.cos(lon) * math.cos(lat)
    y = math.sin(lon) * math.cos(lat)
    z = math.sin(lat)
    return [x,y,z]

def m3x3rotVth(v, theta):
    """Return a 3x3 rotation matrix for a rotation of *theta* radians
    around vector *v*.  In a right-handed systen, positive *theta* gives
    a counter-clockwise rotation when viewed along the line from *v* to
    the origin.
    """

    import math

    # From
    # http://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle

    c = math.cos(theta)
    s = math.sin(theta)
    d = 1-c
    x,y,z = v

    R = [[c+x**2*d,  x*y*d-z*s, x*z*d+y*s],
         [y*x*d+z*s, c+y**2*d,  y*z*d-x*s],
         [z*x*d-y*s, z*y*d+x*s, c+z**2*d]]
    return R

def m3x3multV(M, v):
    """Multiply the vector *v* by the 3x3 matrix *M*: result is Mv."""

    return [VdotV(row,v) for row in M]

def VdotV(u, v):
    """Dot product of *u* and *v*."""

    return sum(x*y for x,y in zip(u,v))


def commands(arg):
    """[command] Show command list."""

    prefix = re.compile(r'\[command] *')

    for name,f in members.items():
        doc = getattr(f, '__doc__')
        if doc and prefix.match(doc):
            doc = prefix.sub('', doc)
            doc = re.sub(r'\n *', ' ', doc)
            doc = re.sub(r'[;.] +.*', '.', doc)
            print name, ':', doc

members = locals()

def main(argv=None):
    import sys
    if argv is None:
        argv = sys.argv

    arg = argv[1:]
    if not arg:
        # With no args, run "commands" to show the command list.
        arg = ['commands']
    if arg[0] in members:
        members[arg[0]](arg)
    else:
        raise Fatal("Unknown command %s." % arg[0])

if __name__ == '__main__':
    main()

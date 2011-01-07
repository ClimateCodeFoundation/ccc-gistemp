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
import re

class Fatal(Exception):
    """An error occurred."""

def cellstations(arg):
    """[command] Compares the lists of stations for each cell.  arg[1]
    and arg[2] should be the names of two different log/step3.log files.
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
            if ina:
                print "Cell %s %s only " % (cell, fs[0].name,
                  ' '.join(sorted(ina)))
            if inb:
                print "Cell %s %s only " % (cell, fs[1].name,
                  ' '.join(sorted(inb)))
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
    from cell identifier (12 characters) to a list of (station,weight)
    pairs.
    """

    result = {}
    for row in f:
        row = row.split(' ', 2)
        if row[1] != 'stations':
            continue
        stations = json.loads(row[2])
        result[row[0]] = stations
    return result


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

    rx = re.compile(r'\[command] *')

    for name,f in members.items():
        doc = getattr(f, '__doc__')
        if doc and rx.match(doc):
            doc = rx.sub('', doc)
            doc = doc.replace('\n', '')
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

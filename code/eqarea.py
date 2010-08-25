#!/usr/bin/env python
# $URL$
# $Rev$
#
# eqarea.py
#
# David Jones, Ravenbrook Limited.

"""Routines for computing an equal area grid.

Specifically, Sergej's equal area 8000 element grid of the Earth (or
any sphere like object).  See GISTEMP code, subroutine GRIDEA

Where "tuple" is used below, some "tuples" may in fact be other
sequence types, such as lists.  Duck typing rules.

"""
__docformat__ = "restructuredtext"

# http://www.python.org/doc/2.3.5/lib/module-itertools.html
import itertools
import math

#: Array of band altitudes.
#:
#: Altitude refers to the distance from the equatorial plane.  The
#: "North" pole being at altitude 1, the "South" pole being at altitude
#: -1.
#:
#: band_altitude[n]
#:   gives the northern altitude of band n, n from 0 to 3.
#:
#: band_altitude[n+1]
#:   gives the southern altitude of band n.
#:
#: Note: These are the sines of the latitude.
band_altitude = [1, 0.9, 0.7, 0.4, 0]

#: Number of horizontal boxes in each band.
#:
#: To ensure equal area these should be proportional to the thickness of each
#: band (Archimedes' Hat Box Theorem).
band_boxes = [4, 8, 12, 16]

# Please move elsewhere
def lerp(x, y, p):
    """Interpolate between x and y by the fraction p.

    When p == 0 x will be returned, when p == 1 y will be returned.  Note
    that p is not restricted to being between 0 and 1.

    :Return:
        The interpolated value.

    :Param x, y:
        The interpolation end-points.
    :Param p:
        The interpolation fraction.

    """
    p = float(p)
    return y*p + (1-p)*x


def northern40() :
    """Generator: Yields the 40 northern hemisphere boxes.

    The yielded value is a tuple of::

        (southern, northern, western, eastern).

    See `grid` for more details about the coordinates.

    """
    for band in range(len(band_boxes)) :
        # number of horizontal boxes in band
        n = band_boxes[band]
        for i in range(n) :
            lats = 180/math.pi*math.asin(band_altitude[band+1])
            latn = 180/math.pi*math.asin(band_altitude[band])
            lonw = -180 + 360*float(i)/n
            lone = -180 + 360*float(i+1)/n
            yield (lats, latn, lonw, lone)


def southern40() :
    """Generator: Yields the 40 southern hemisphere boxes.i

    The yielded value is a tuple of::

        (southern, northern, western, eastern).

    See `grid` for more details about the coordinates.

    """
    # Note: Avoid "reversed" because it is not in Python 2.3
    n = list(northern40())
    # We want to take the northern list and reverse the bands, within each
    # band the boxes will still be in the same order.  So first
    # gather into bands.
    i = 0
    band = []
    for w in band_boxes:
        band.append(n[i:i+w])
        i += w
    # Then reverse the band list
    band.reverse()
    # And stitch back into a single, southern, list
    s = []
    # Note: used for side-effect!
    map(lambda x: s.extend(x), band)
    assert len(s) == len(n)

    # Flip each box north/south
    for x in s :
        yield (-x[1],-x[0],x[2],x[3])


def grid() :
    """Generator: Yields a list of 80 boxes equally dividing a sphere.

    Each box comprises an equal area division of a sphere; each box is
    described by a 4-tuple of its boundaries::

        (southern, northern, western, eastern).

    Co-ordinates are given as (fractional) degrees of latitude (for
    northern and southern borders) and longitude (for western and eastern
    borders).

    """
    return itertools.chain(northern40(), southern40())


def gridsub() :
    """Generator: Yields 80 boxes each containing a subbox generator.

    Each yielded box contains, in a general sense, 100 subboxes.
    The 80 boxes are those corresponding to grid().  A box is returned as a
    pair (bound, subgen), where bound is the 4-tuple as returned by `grid`
    (boundaries in degrees for S N W E edges), and subgen is a generator for
    the 100 subboxes.

    The order of the boxes is the same as `grid`.

    The subboxes are returned as 4-tuples using the same
    latitude/longitude box representation as `grid`.

    Note that Sheffield, +53.40-001.50, is in subbox 759:

    >>> subbox=list(list(gridsub())[7][1])[59]
    >>> subbox[0] < 53.4 < subbox[1]
    True
    >>> subbox[2] < -1.5 < subbox[3]
    True

    and Vostok, -78.40+106.90, is in subbox 7921:

    >>> subbox=list(list(gridsub())[79][1])[21]
    >>> subbox[0] < -78.4 < subbox[1]
    True
    >>> subbox[2] < 106.9 < subbox[3]
    True

    """
    def subgen(box):
        """A generator for the subboxes of box."""

        # Altitude for southern and northern border.
        alts = math.sin(box[0]*math.pi/180)
        altn = math.sin(box[1]*math.pi/180)
        for y in range(10):
            s = 180*math.asin(lerp(alts, altn, y*0.1))/math.pi
            n = 180*math.asin(lerp(alts, altn, (y+1)*0.1))/math.pi
            for x in range(10):
                w = lerp(box[2], box[3], x*0.1)
                e = lerp(box[2], box[3], (x+1)*0.1)
                yield(s, n, w, e)

    for box in grid():
        yield (box, subgen(box))


def grid8k() :
    """Generator: As `gridsub`, but flattened.

    Yields the same set of boxes as `gridsub`, but returns a single generator
    for all 8000 subboxes.  Not used by core code, but used by tools and useful
    for debugging.

    """
    for box in gridsub():
        for subbox in box[1]:
            yield subbox


def gridR3() :
    """Generator: Yields a list of 80 R3-coordinate boxes equally dividing a
    sphere.

    Like `grid` but each box is described in a right-hand 3-dimensional
    euclidean co-ordinate system.  The sphere is centred at the origin
    ``(0,0,0)``. The equator lies on the plane containing the x- and y-axes
    ``(z=0)``.  The north pole has co-ordinate ``(0,0,1)``; the intersection of
    the eqautor and the prime meridian (0 longitude at equator) is
    ``(1,0,0)``.  That should be enough for you to deduce that longitude 90 on
    the equator is at ``(0,1,0)``.  Each box is a quadrilateral described by
    its corners; from the viewpoint of a distant observer viewing the side
    of the sphere with the box on, the corners are in a counter-clockwise
    order (so this consistent ordering can be used for face culling if
    desired).  Each box is described by a 4-tuple of 3-tuples.  Note polar
    boxes have co-incident corners.

    """
    def llto3d(i) :
      """Convert (latitutude, longitude) into 3D triple (x,y,z).  Latitude
      and longitude are given in degrees."""

      lat = i[0]*math.pi/180
      z = math.sin(lat)
      c = math.cos(lat)
      long = i[1]*math.pi/180
      x = math.cos(long) * c
      y = math.sin(long) * c
      return (x,y,z)

    def bto3d(i) :
      """Convert 4-tuple of borders into counter-clockwise 4-tuple of 3d
      co-ordinates."""

      # For counter-clockwise corners, start at NE and work round.
      # Recall border 4-tuple convention described in grid.
      return map(llto3d, ((i[1],i[3]), (i[1],i[2]), (i[0],i[2]), (i[0],i[3])))

    return itertools.imap(bto3d, grid())


def gridJSON() :
    """Create array of 80 R3-coordinate boxes in JSON format.

    As `gridR3` but in JSON (http://www.json.org/) format.  Returned
    string matches a JSON array production.

    """
    return str(list(itertools.imap(lambda x: map(list, x), gridR3())))


def centre(box) :
    """Calculate the (latitude,longitude) pair for the centre of box/subbox.

    This the the "equal area" centre in the sense that the area to the
    north will equal the area to the south, and the same for east/west.

    :Return:
        The ``(latitude,longitude)`` for the box or subbox.

    :Param box:
        The box (or subbox) to find the centre of. Specified as a 4-tuple
        of its boundaries (same convention used by grid()):
        (southern, northern, western, eastern).

    """
    sinc = 0.5*(math.sin(box[0]*math.pi/180) + math.sin(box[1]*math.pi/180))
    return (math.asin(sinc)*180/math.pi, 0.5*(box[2]+box[3]))


def main() :
    # http://www.python.org/doc/2.3.5/lib/module-doctest.html
    import doctest, eqarea
    return doctest.testmod(eqarea)


if __name__ == '__main__' :
    main()

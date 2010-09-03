#!/usr/bin/env python
# $URL$
# $Rev$

"""
Convert GHCN style data and metadata file to a KML file for browsing
in Google Earth.

Usage:
python tool/kml.py [--longitude west,east] [--latitude south,north]
  [data [inv]] [> ghcn.kml]

The data file is usually v2.mean, it defaults to input/v2.mean; the
inv file, metadata, defaults to input/v2.inv.

Using --longitude or --latitude restricts the output to stations within
those bounds.
"""

# Clear Climate Code
import extend_path
# Clear Climate Code
from giss_io import MISSING


# Images for the icons that denote the "urban flag" in the v2.mean file;
# this flag is either R (rural), S (semi-urban), U (urban).
urbanicon = dict(zip('RSU', """
http://google-maps-icons.googlecode.com/files/forest.png
http://google-maps-icons.googlecode.com/files/home.png
http://google-maps-icons.googlecode.com/files/factory.png
""".split()))

def xmlquote(x):
    import re

    def it(m):
        """Takes match object, returns replacement string."""
        s = m.group()
        if s == '&':
            return '&amp;'
        if s == '<':
            return '&lt;'

    return re.sub('[&<]', it, x)

def meta_iter(meta):
    """Convert v2.mean metadata, in a v2.inv file, into a metadata
    iterator."""

    for line in meta:
        yield line[:11], line

# Regular expression for a record identifer for a cell.  These are ISO
# 6709 lat/lon value in degrees with a precision of .1.
CELL_RE = r'([-+]\d\d\.\d)([-+]\d\d\d\.\d)'


def kml(data, meta, out, lat=(-90,91), lon=(-180,181)):
    """Convert the input *meta*, which is assumed to be (an open file
    object on) the NOAA file v2.temperature.inv file or one in an
    equivalent format, to a KML file.  The output is written to *out*.

    The KML file can be opened with Google Earth.
    """

    import re

    # Make meta a dict as we need to use it several times.
    meta = dict(meta_iter(meta))

    # Most of the KML knowledge was scraped from
    # http://code.google.com/apis/kml/documentation/kml_tut.html
    # and
    # http://code.google.com/apis/kml/documentation/kmlreference.html
    print '<?xml version="1.0" encoding="UTF-8"?>'
    print '<kml xmlns="http://www.opengis.net/kml/2.2">'
    print '<Document>'

    annual_data = station_annual_anomalies(data)
    styled = set()
    for id, anom in annual_data:
        # Bit hacky, just use first record when multiple records.
        id11 = id[:11]
        if id11 not in styled:
            m = re.match(CELL_RE, id11)
            if m:
                # Grid cell.
                clat,clon = m.groups()
                # Create a fake metadata string
                # Calculate first year with data (with yearly anomaly
                # really).
                first = [i for i,v in enumerate(anom) if v != MISSING]
                if first:
                    first = str(first[0]+1880)
                else:
                    first = ''
                name = first
                fake = "%(id11)s %(name)-30s %(clat)-6s %(clon)-7s" % locals()
                fake = "%-80s" % fake
                meta[id11] = fake
            # Assume normal 11-digit GHCN station.
            # todo, use metadata here, for urban/airport etc.
            url = chart_anom(anom)
            print """
<Style id="station%(id)s">
  <IconStyle><Icon>
    <href>%(url)s</href>
  </Icon></IconStyle>
</Style>""" % dict(id=id11, url=xmlquote(url))
            styled.add(id11)

    for stationid,l in meta.items():
        # See ftp.ncdc.noaa.gov/v2.read.inv.f
        # A 30 character name field.  It often contains a country and
        # lots of internal spaces.  We reduce all multiples of SPACE
        # character to two SPACE characters.
        name = l[12:42]
        name = name.strip()
        name = re.sub(' {2,}', '  ', name)
        name = xmlquote(name)

        # v2.inv has north then east, KML has east then north.
        north = l[43:49]
        east = l[50:57]
        north = north.strip()
        east = east.strip()
        if not (lat[0] <= float(north) < lat[1] and
          lon[0] <= float(east) < lon[1]):
            continue
        # (no longer) used to determine icon (see urbanicon, above)
        urbanflag = l[67]

        if stationid not in styled:
            continue

        print '<Placemark>'
        print '<name>%(name)s  %(stationid)s</name>' % locals()
        print '<description>%s</description>' % stationid
        print '<Point><coordinates>'
        print "%(east)s,%(north)s" % locals()
        print '</coordinates></Point>'
        print '<styleUrl>#station%s</styleUrl>' % stationid
        print '</Placemark>'
    print '</Document>'
    print '</kml>'

def station_annual_anomalies(data):
    """*data* should be a v2.mean file.  Return a dictionary that
    maps from station record identifier (12-digit) to an annual anomaly
    series.  All the data series start in 1880."""

    from code import step1
    from giss_io import V2MeanReader

    return ((record.uid, step1.monthly_annual(record.series)[1])
      for record in V2MeanReader(data, year_min=1880))

def chart_anom(anom):
    """Generate a (Google Chart API) URL for a chart of annual
    anomalies."""

    prefix = "http://chart.apis.google.com/chart"

    def treat(datum):
        """Convert datum to string form."""

        if datum == MISSING:
            return '-999'
        return '%.2f' % datum

    series = ','.join(map(treat, anom))

    chart = [
      'chs=64x64',
      'chd=t:' + series,
      'chds=-5,5',
      'cht=ls',
      'chco=ffffff',
      'chf=bg,s,FF0000',
    ]

    return prefix + '?' + '&'.join(chart)

def usage():
    import sys
    print >>sys.stderr, __doc__

def main(argv=None):
    import sys
    import getopt

    if argv is None:
        argv = sys.argv

    key = {}
    try:
        opts, arg = getopt.getopt(argv[1:], '', ['longitude=', 'latitude='])
    except getopt.GetoptError:
        return usage()
    for o,v in opts:
        if o == '--longitude':
            key['lon'] = map(float, v.split(','))
        elif o == '--latitude':
            key['lat'] = map(float, v.split(','))

    if arg:
        d = arg[0]
    else:
        d = 'input/v2.mean'
    if len(arg) > 1:
        m = open(arg[1])
    else:
        m = open('input/v2.inv')
    return kml(d, m, sys.stdout, **key)

if __name__ == '__main__':
    main()

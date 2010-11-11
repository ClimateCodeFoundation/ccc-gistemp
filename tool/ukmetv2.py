#!/usr/bin/env python
# $URL$
# $Rev$
#
# ukmetv2.py
#
# Scrape and convert UK Met Office station data from their website and
# convert to GHCN v2.mean format.
#
# Requires: Python 2.5 (or higher) and BeautifulSoup.

import itertools
# http://docs.python.org/release/2.5.4/lib/module-math.html
import math
import re
import sys
# http://docs.python.org/release/2.4.4/lib/module-traceback.html
import traceback
import urllib

# Requires Python 2.5.
assert sys.version_info[:2] >= (2,5)

# Clear Climate Code
import extend_path
from code.giss_data import MISSING
import io

PAGE = """http://www.metoffice.gov.uk/climate/uk/stationdata/"""

class Struct:
    pass

# Static list of pages, instead of scraping index page.  Used when
# BeautifulSoup is not available.
static_list = """
http://www.metoffice.gov.uk/climate/uk/stationdata/aberporthdata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/armaghdata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/ballypatrickdata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/bradforddata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/braemardata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/cambornedata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/cambridgedata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/cardiffdata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/chivenordata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/cwmystwythdata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/dunstaffnagedata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/durhamdata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/eastbournedata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/eskdalemuirdata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/heathrowdata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/hurndata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/lerwickdata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/leucharsdata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/lowestoftdata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/manstondata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/nairndata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/newtonriggdata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/oxforddata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/paisleydata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/ringwaydata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/rossonwyedata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/shawburydata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/sheffielddata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/southamptondata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/stornowaydata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/suttonboningtondata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/tireedata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/valleydata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/waddingtondata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/whitbydata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/wickairportdata.txt
http://www.metoffice.gov.uk/climate/uk/stationdata/yeoviltondata.txt
""".split()

def scrapeit(prefix):
    """Creates files of the form:
      prefix.v2.mean
      prefix.v2.min
      prefix.v2.max
      prefix.v2.inv
    """

    out = Struct()
    out.mean = io.GHCNV2Writer(path=prefix + '.v2.mean')
    out.min = io.GHCNV2Writer(path=prefix + '.v2.min')
    out.max = io.GHCNV2Writer(path=prefix + '.v2.max')
    out.inv = open(prefix + '.v2.inv', 'w')

    for url in geturls():
        print "from", url
        scrape1(url, out)

def geturls():
    """Return a sequence of URLs, each one being a text file containing
    temperature data.  Scraped from "Historic station data" of the UK
    Met Office."""

    # http://docs.python.org/release/2.4.4/lib/module-urlparse.html
    import urlparse

    try:
        from BeautifulSoup import BeautifulSoup
    except ImportError:
        BeautifulSoup = None

    if not BeautifulSoup:
        for u in static_list:
            yield u
        return

    text = urllib.urlopen(PAGE).read()
    souped = BeautifulSoup(text)

    for element in souped.findAll('option', value=re.compile(r'\.txt$')):
        u = element['value']
        url = urlparse.urljoin(PAGE, u)
        yield url

LARGE = 8888

MONTHS = """January February March April May June July August September
October November December""".split()
MONTHRE = '|'.join(["(%s[a-z]*)" % mon[:3] for mon in MONTHS])
MONTHRE = '(?:' + MONTHRE + ')'

def scrape1(url, out):
    """Scrape one url from the met office website.  URL contains
    temperature data, for example:
    http://www.metoffice.gov.uk/climate/uk/stationdata/yeoviltondata.txt
    """

    gotmeta = False
    f = urllib.urlopen(url)
    name = f.readline()
    try:
        nl = name.split('/')
        if len(nl) > 1:
            print "Split location"
        shortname = re.sub(r'\s+', '', name)
        shortname = (shortname + '_'*6)[:7]
        assert 7 == len(shortname)
        id11 = "UKMO" + shortname.upper()
        id12 = id11 + '0'
        location = f.readline()
        assert location.lower().startswith('location')

        if 0:
            fromre = (
              r'from\s*' + MONTHRE + '\s*(\d{4})')
            fromre = re.compile(fromre, re.IGNORECASE)
            froms = re.findall(fromre, location)
            if froms:
                print "froms", froms

        # Parse out the grid location (or locations, some stations have
        # moved).
        # We find each section of string that is a grid position, then
        # we try and parse the bits following a grid position for a date
        # range.
        # Grid position is generally of the form:
        # 4420E 1125N, 20 metres amsl
        # But optional commas, spaces, "m"/"metres".  And two stations are
        # on the island of Ireland and therefore use the Irish Grid.
        GRIDRE = re.compile(
          r'(\d+)E\s*(\d+)N(\s*\(Irish Grid\))?(?:\s|,)*(\d+)\s*m')
        gridmatches = list(re.finditer(GRIDRE, location))
        # In *locs* build a list of objects that have a .grid member for
        # the match object, and a .rest member for the remainder of the
        # string up to the next match.
        locs = []
        for i,m in enumerate(gridmatches):
            loc = Struct()
            # The match object for GRIDRE
            loc.grid = m
            if i == len(gridmatches)-1:
                loc.rest = location[m.end():]
            else:
                loc.rest = location[m.end():gridmatches[i+1].start()]
            locs.append(loc)
        if len(locs) != 1:
            print len(locs), "grid locations found"
        meta = Struct()
        old = None
        for loc in locs:
            grid = loc.grid.groups()
            easting, northing, irishgrid, height = grid
            if len(locs) > 1:
                print ','.join([easting, northing, height]), grid[2]
            # Convert to number, metres.
            easting, northing, height = map(int, [easting, northing, height])
            # Easting and Northing given are in units of 100m.  It
            # doesn't say that, but they are.  Note in particular
            # http://www.metoffice.gov.uk/climate/uk/stationdata/cambornedata.txt
            easting, northing = [100*x for x in [easting, northing]]
            if old:
                print "Move: %.0f along; height change: %+d" % (
                  math.hypot(easting-old[0], northing-old[1]),
                  height-old[2])
            old = (easting, northing, height)
            # Don't forget to use http://gridconvert.appspot.com/ to convert
            # from OSGB to GRS80. (which won't work for Irish Grid)
            irishgrid = grid[2]
            if irishgrid:
                print "Irish Grid, can't compute lat/lon."
            else:
                lat,lon,h = osgbtolatlon(easting, northing, height)
                meta.lat = lat
                meta.lon = lon
        meta.uid = id12
        meta.name = nl[0].strip()
        # Note: Height above Mean Sea Level, not GRS80 ellipsoid height.
        meta.height = height
        out.inv.write(metav2(meta) + '\n')

    except AssertionError, err:
        print "Skipping metadata"
        traceback.print_exc()

    def year(row):
        """Extract the year from a row."""
        return re.match(r'\s*(\d{4})', row).group(1)

    # List of splits created by hand
    # Key is (uid, year, month) triple, value is the new 12 character
    # uid (which sticks until there is another split).  year, month
    # specifies the first month of the newly split record.
    split = {
      ('UKMOWHITBYC0',2000,1): 'UKMOWHITBY_$',
      }

    for y,rows in itertools.groupby(datarows(f), year):
        # Split at beginning of year.
        yr = int(y)
        newsplit = split.get((id12, yr, 1), None)
        if newsplit:
            id12 = newsplit
            print "Changing ID to", id12
        # 12 months of temps:
        mins = [LARGE]*12
        maxs = [LARGE]*12
        means = [LARGE]*12
        for row in rows:
            m,tmin,tmax = row.split()[1:4]
            m = int(m) - 1
            mins[m] = convert1(tmin)
            maxs[m] = convert1(tmax)
        means = [(tmin+tmax)/2.0 for tmin,tmax in zip(mins,maxs)]
        means = [x if x < 1000 else MISSING for x in means]
        mins = [x if x < 1000 else MISSING for x in mins]
        maxs = [x if x < 1000 else MISSING for x in maxs]
        out.mean.writeyear(id12, yr, means)
        out.min.writeyear(id12, yr, mins)
        out.max.writeyear(id12, yr, maxs)

def osgbtolatlon(easting, northing, height):
    """*easting* and *northing* are the easting and northing in metres
    on the Ordanance Survery National Grid.  *height* is the height
    above mean sea level.  Return (lat,lon,height) triple for the GRS80
    ellipsoid (returned height is ellipsoid height).
    """

    # See http://gridconvert.appspot.com/
    SERVICE = "http://gridconvert.appspot.com/osgb36/etrs89/"
    url = "%s%r,%r,%r" % (SERVICE, easting, northing, height)
    result = urllib.urlopen(url).read()
    # *result* is a JSON object, but we parse it with regular
    # expressions (to avoid a json dependency).
    # RE for a JSON number.  See http://www.json.org/
    json_no = r'(-?\d+(?:\.\d+)?)'
    try:
        lat = re.search(r'ETRS89_Latitude[^-\d]*' + json_no, result).group(1)
        lon = re.search(r'ETRS89_Longitude[^-\d]*' + json_no, result).group(1)
        height = re.search(r'ETRS89_Height[^-\d]*' + json_no, result).group(1)
    except:
        return None
    return float(lat),float(lon),height

def metav2(meta):
    """Convert a metadata object for a station into a v2.inv style
    record.  A string is returned (with no trailing newline).
    """

    # Perhaps the best documentation for the format is the `stations`
    # function in code/giss_data.py

    inv = [' ']*106
    def convertlat(x):
        return "%6.2f" % x
    def convertlon(x):
        return "%7.2f" % x
    def convertheight(x):
        return "%4d" % x
    fields = (
        (0,   11,  'uid',               str),
        (12,  42,  'name',              str),
        (43,  49,  'lat',               convertlat),
        (50,  57,  'lon',               convertlon),
        (58,  62,  'height',            convertheight),
    )
    for a,b,attr,convert in fields:
        if hasattr(meta, attr):
            item = getattr(meta, attr)
            item = convert(item)
            # Pad/truncate to length.
            item = item.ljust(b-a)[:b-a]
            inv[a:b] = item
    return ''.join(inv)

def convert1(t):
    """Convert 1 temperature datum from a scraped file.  These are
    decimal temperature values (in degrees C).  Absent data are marked
    with "---"; estimated data are followed by a "*".  Currently,
    estimated data are discarded.  If the datum is absent (or
    estimated), LARGE is returned.
    """

    if t.endswith('*'):
        return LARGE
    if '---' in t:
        return LARGE
    return float(t)
    

def datarows(f):
    """Filter a stream so that only rows of data are yielded.  This
    yields every row that starts with a 4 digit year (and possibly some
    whitespace)."""

    for row in f:
        if re.match(r'\s*\d{4}', row):
            yield row
    

def main(argv=None):
    import sys
    if argv is None:
        argv = sys.argv

    scrapeit("ukmet")

if __name__ == '__main__':
    main()

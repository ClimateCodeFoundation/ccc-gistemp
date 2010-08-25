#!/usr/bin/env python
# $URL$
# $Rev$
#
# ukmetv2.py
#
# Scrape and convert UK Met Office station data from their website and
# convert to GHCN v2.mean format.

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
import giss_io

PAGE = """http://www.metoffice.gov.uk/climate/uk/stationdata/"""

class Struct:
    pass

def scrapeit(prefix):
    """Creates files of the form:
      prefix.v2.mean
      prefix.v2.min
      prefix.v2.max
      prefix.v2.inv
    """

    out = Struct()
    out.mean = giss_io.V2MeanWriter(path=prefix + '.v2.mean')
    out.min = giss_io.V2MeanWriter(path=prefix + '.v2.min')
    out.max = giss_io.V2MeanWriter(path=prefix + '.v2.max')
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

    from BeautifulSoup import BeautifulSoup

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
        old = None
        for loc in locs:
            grid = loc.grid.groups()
            # Extend Easting and Northing to 6 figures (metres)
            easting, northing = [x + '0'*(6-len(x)) for x in grid[0:2]]
            height = grid[3]
            if len(locs) > 1:
                print ','.join([easting, northing, height]), grid[2]
            # Convert to number, metres.
            easting, northing, height = map(int, [easting, northing, height])
            if old:
                print "Move: %.0f along; height change: %+d" % (
                  math.hypot(easting-old[0], northing-old[1]),
                  height-old[2])
            old = (easting, northing, height)
            # Don't forget to use http://gridconvert.appspot.com/ to convert
            # from OSGB to GRS80. (which won't work for Irish Grid)
            # :todo: save metadata instead of ignoring it.
        meta = Struct()
        meta.uid = id12
        meta.name = nl[0].strip()
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

def metav2(meta):
    """Convert a metadata object for a station into a v2.inv style
    record.  A string is returned (with no trailing newline).
    """

    # Perhaps the best documentation for the format is the `stations`
    # function in code/giss_data.py

    inv = [' ']*106
    fields = (
        (0,   11,  'uid',               str),
        (12,  42,  'name',              str),
        (58,  62,  'height',            str),
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

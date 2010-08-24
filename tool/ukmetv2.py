#!/usr/bin/env python
# $URL$
# $Rev$
#
# ukmetv2.py
#
# Scrape and convert UK Met Office station data from their website and
# convert to GHCN v2.mean format.

import itertools
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

def scrapeit(prefix):
    """Creates files of the form:
      prefix.v2.mean
      prefix.v2.min
      prefix.v2.max
      prefix.v2.inv
    """

    class Struct:
        pass

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

def scrape1(url, out):
    """Scrape one url from the met office website.  URL contains
    temperature data, for example:
    http://www.metoffice.gov.uk/climate/uk/stationdata/yeoviltondata.txt
    """

    meta = False
    f = urllib.urlopen(url)
    name = f.readline()
    try:
        nl = name.split('/')
        if len(nl) > 1:
            print "Split location"
        shortname = name.replace(' ', '')
        shortname = (shortname + '_'*6)[:7]
        assert 7 == len(shortname)
        id11 = "UKMO" + shortname.upper()
        id12 = id11 + '0'
        location = f.readline()
        assert location.lower().startswith('location')
        # Optional commas, spaces, "m"/"metres".  And two stations are
        # on the island of Ireland and therefore use the Irish Grid.
        GRIDRE = re.compile(
          r'(\d+)E\s*(\d+)N(\s*\(Irish Grid\))?(?:\s|,)*(\d+)\s*m')
        locs = re.findall(GRIDRE, location)
        if len(locs) != 1:
            print len(locs), "grid locations found"
        for loc in locs:
            # Extend Easting and Northing to 6 figures (metres)
            loc = [x + '0'*(6-len(x)) for x in loc[0:2]] + list(loc[2:])
            # Convert to number, metres.
            # easting, northing, irishgrid, height = map(int, loc)
            # Don't forget to use http://gridconvert.appspot.com/ to convert
            # from OSGB to GRS80. (which won't work for Irish Grid)
            # :todo: save metadata instead of ignoring it.
            if len(locs) > 1:
                print ','.join(loc)

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
      ('UKMOWHITBYC0',2000,1): 'UKMETOWHITBY_',
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

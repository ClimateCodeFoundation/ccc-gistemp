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
import urllib

# Requires Python 2.5.
assert sys.version_info[:2] >= (2,5)

# Clear Climate Code
import extend_path
from code.giss_data import MISSING
import giss_io

urls="""
http://www.metoffice.gov.uk/climate/uk/stationdata/oxforddata.txt
""".split()

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


    for url in getlocations():
        print "Data from", url, "..."
        scrape1(url, out)

def getlocations():
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
        shortname = name.replace(' ', '')
        shortname = (shortname + '_'*6)[:7]
        assert 7 == len(shortname)
        id11 = "UKMO" + shortname.upper()
        id12 = id11 + '0'
        location = f.readline()
        assert location.lower().startswith('location')
        loc = location.split()[1:3]
        loc[1] = loc[1].strip(',')
        # Sort by last char, so that easting comes first.
        loc.sort(key=lambda x:x[-1])
        assert loc[0].endswith('E')
        assert loc[1].endswith('N')
        loc = [x[:-1] for x in loc]
        # Extend to 6 figures (metres)
        loc = [x + '0'*(6-len(x)) for x in loc]
        # Convert to number, metres.
        easting, northing = map(int, loc)
        # Don't forget to use http://gridconvert.appspot.com/ to convert
        # from OSGB to GRS80.
        # :todo: save metadata instead of ignoring it.
    except AssertionError, err:
        print "Skipping metadata"
        print err

    def year(row):
        """Extract the year from a row."""
        return re.match(r'\s*(\d{4})', row).group(1)

    for y,rows in itertools.groupby(datarows(f), year):
        yr = int(y)
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

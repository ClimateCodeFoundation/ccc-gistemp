#!/usr/bin/env python
# $URL$
# $Rev$

# cascrape.py
#
# David Jones, Climate Code Foundation, 2010-10-26

"""
Tool to collect temperature data for Candada scraped by ScraperWiki.
"""

import json
import urllib

def gather(scraper, v2out, v2inv, log):
    """Takes the name of a scraper hosted at ScraperWiki and gathers all
    its data."""

    api = "http://api.scraperwiki.com/api/1.0/datastore/getdata"
    name = scraper

    # We are limited by the API to retrieving a certain number of
    # records per request, so we loop through, increasing offset, until
    # we have got them all.
    limit = 500
    offset = 0
    while 1:
        url = api + '?' + urllib.urlencode(dict(name=name,
            format='json', limit=limit, offset=offset))
        gotnone = True
        result = json.load(urllib.urlopen(url))
        for item in result:
            if 'meta' in item:
                log.write(item['meta'])
                log.write('\n')
                meta = json.loads(item['meta'])
                v2inv.write(formatmetav2(meta))
                continue
            # Make a row in GHCN v2 format.  Add a '0' to make a
            # 12-character record identifier.
            row = "%s0%s%s\n" % (item['id'], item['year'], item['tmean'])
            v2out.write(row)
            gotnone = False
        if gotnone:
            break
        offset += 500
        # Just a sanity check
        assert offset < 1e5

def formatmetav2(meta):
    """Format a metadata dict as a row in the same format as the GHCN
    v2.temperature.inv file."""

    assert 11 == len(meta['id11'])
    for key in ['Latitude', 'Longitude', 'Elevation']:
        meta[key] = float(meta[key])
    invrow = (
      "%(id11)s "
      "%(Station Name)-30s "
      "%(Latitude)6.2f "
      "%(Longitude)7.2f "
      "%(Elevation)4.0f "
      "-999            -9 -9         \n"
      ) % meta
    return invrow


def main(argv=None):
    import sys
    if argv is None:
        argv = sys.argv

    gather('canada-climate-data',
      v2out=open('ca.v2', 'w'),
      v2inv=open('ca.v2.inv', 'w'),
      log=sys.stdout)

if __name__ == '__main__':
    main()

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

def gather(scraper, v2out, v2inv, jsoninv, log):
    """Takes the name of a scraper hosted at ScraperWiki and gathers all
    its data."""

    api = "http://api.scraperwiki.com/api/1.0/datastore/getdata"
    name = scraper

    # We are limited by the API to retrieving a certain number of
    # records per request, so we loop through, increasing offset, until
    # we have got them all.
    limit = 500
    offset = 0
    meta = []
    rows = []
    count = 0
    while 1:
        url = api + '?' + urllib.urlencode(dict(name=name,
            format='json', limit=limit, offset=offset))
        # We terminate when the JSON request returns no results, which
        # we detect using this flag.
        gotnone = True
        result = json.load(urllib.urlopen(url))
        for item in result:
            gotnone = False
            count += 1
            log.write('\r%d records' % count)
            log.flush()
            if 'meta' in item:
                meta.append(json.loads(item['meta']))
                continue
            # Make a row in GHCN v2 format.  Add a '0' to make a
            # 12-character record identifier.
            id11 = item['id']
            assert 11 == len(id11)
            row = "%s%s%s\n" % (id11+'0', item['year'], item['tmean'])
            rows.append(row)
        if gotnone:
            break
        offset += 500
        # Just a sanity check
        assert offset < 1e5
    log.write('\n')
    v2out.writelines(sorted(rows))
    for m in sorted(meta, key=lambda x: x['id11']):
        jsoninv.write(json.dumps(m) + '\n')
        v2inv.write(formatmetav2(m))

def formatmetav2(meta):
    """Format a metadata dict as a row in the same format as the GHCN
    v2.temperature.inv file."""

    assert 11 == len(meta['id11'])
    for key in ['Latitude', 'Longitude', 'Elevation']:
        meta[key] = float(meta[key])
    invrow = (
      "%(id11)s "
      "%(Station Name)-30.30s "
      "%(Latitude)6.2f "
      "%(Longitude)7.2f "
      "%(Elevation)4.0f "
      "-999            -9 -9UNKNOWN         "
      "      "
      "\n"
      ) % meta
    return invrow


def main(argv=None):
    import sys
    if argv is None:
        argv = sys.argv

    gather('canada-climate-data',
      v2out=open('input/ca.v2.mean', 'w'),
      v2inv=open('input/ca.v2.inv', 'w'),
      jsoninv=open('input/ca.json', 'w'),
      log=sys.stdout)

if __name__ == '__main__':
    main()

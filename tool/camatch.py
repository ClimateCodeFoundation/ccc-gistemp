#!/usr/bin/env python
# $URL$
# $Rev$
#
# camatch.py
#
# Match stations in the scraped Canada data (see cascrape.py) with GHCN
# stations.

import json

# Clear Climate Code
import extend_path
from code import giss_data

def match(cameta, ghcnmeta):
    """Analyse and match scraped Canada stations with GHCN stations."""

    cadict = cametaasdict(cameta)
    # 403 is the GHCN prefix for Canada, see v2.country.codes
    ghcndict = ghcnmetaasdict(ghcnmeta, prefix='403')

    wmocount = 0
    wmomatch = 0
    wmolocnear = 0
    wmolocfail = 0
    locmatch = 0
    nomatch = 0
    for cast in cadict.values():
        wmo = cast['WMO Identifier']
        caloc = iso6709(cast['Latitude'], cast['Longitude'])
        if wmo:
            wmocount += 1
            wmoid11 = "403%s000" % wmo
            if wmoid11 in ghcndict:
                wmomatch += 1
                print "match WMO %s %s" % (cast['id11'], wmoid11)
                ghcnst = ghcndict[wmoid11]
                ghcnloc = iso6709(ghcnst.lat, ghcnst.lon)
                if caloc != ghcnloc:
                    if locnear((cast['Latitude'], cast['Longitude']),
                      (ghcnst.lat, ghcnst.lon)):
                          wmolocnear += 1
                    else:
                        wmolocfail += 1
                        print ("location different: WMO %s CA %s GHCN %s" %
                          (wmo, caloc, ghcnloc))
                continue
        if caloc in ghcndict:
            locmatch += 1
            print "match LOCATION %s %s %s" % (
              cast['id11'], ghcndict[caloc].uid, wmo)
            continue
        nomatch += 1

    print "WMO stations", wmocount
    print "WMO match", wmomatch, "(of those) near:", wmolocnear, "not near: ", wmolocfail
    print "LOCATION match", locmatch
    print "no match", nomatch


def cametaasdict(meta):
    """Take an open file descriptor on the 'ca.json' file and return a
    dictionary."""

    def itermeta():
        for row in meta:
            item = json.loads(row)
            for key in ['Latitude', 'Longitude', 'Elevation']:
                item[key] = float(item[key])
            yield item['Climate Identifier'], item
    return dict(itermeta())

def ghcnmetaasdict(ghcnmeta, prefix=None):
    """Take an open file descriptor on the 'v2.inv' file and return a
    dictionary."""

    stations = giss_data.read_stations(file=ghcnmeta, format='v2')
    if prefix:
        keeps = dict((id11,station)
          for id11,station in stations.items()
          if id11.startswith(prefix))
        stations = keeps
    # Additionally index stations by location.
    new = ((iso6709(station.lat, station.lon),station)
      for station in stations.values())
    stations.update(new)
    return stations

def iso6709(lat, lon):
    """Return the lat,lon pair as a string in ISO 6709 format."""

    return "%+06.2f%+07.2f" % (lat, lon)

def locnear(a, b):
    """True when the two points (each a (lat,lon) pair are "near" each
    other.  In this case near means both lat and lon are within 0.01."""
    # Ignores issues at -180/+180
    maxnorm = max(map(lambda x,y: x-y, a, b), key=abs)
    return maxnorm < 0.013

def main(argv=None):
    import sys
    if argv is None:
        argv = sys.argv
    match(cameta=open('ca.json'), ghcnmeta=open('input/v2.inv'))

if __name__ == '__main__':
    main()

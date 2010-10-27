#!/usr/bin/env python
# $URL$
# $Rev$
#
# camatch.py
#
# Match stations in the scraped Canada data (see cascrape.py) with GHCN
# stations.

import itertools
import json
import math

# Clear Climate Code
import extend_path
from code import giss_data

class Struct:
    pass

wmocount = 0

def match(cameta, cadata, ghcnmeta, ghcndata):
    """Analyse and match scraped Canada stations with GHCN stations."""

    cadict = cametaasdict(cameta)
    # 403 is the GHCN prefix for Canada, see v2.country.codes
    ghcndict = ghcnmetaasdict(ghcnmeta, prefix='403')

    drop_short(cadict, cadata)

    global wmocount
    wmomatch = 0
    wmofar = 0
    locnear = 0

    for match in itermatches(cadict, ghcndict):
        cast = match.cast
        ghcnst = match.ghcnst
        sep = locdist((cast['Latitude'], cast['Longitude']),
                      (ghcnst.lat, ghcnst.lon))
        if match.type == 'wmo':
            wmomatch += 1
            if sep > 0.015:
                wmofar += 1
        else:
            if sep < 0.015:
                locnear += 1
        print match.type, match.ghcnst.uid, match.cast['id11'], \
          "%.2f" % sep, match.cast['WMO Identifier']

    print "dropped stations", dropcount
    print "kept stations", len(cadict)
    print "WMO stations", wmocount
    print "WMO match", wmomatch, "(of those) not near:", wmofar
    print "LOCATION near", locnear

def itermatches(cadict, ghcndict):
    """Iterate over all the stations in cadict, yielding a match object
    for each one.  match.type is either 'wmo' or 'near'; match.cast
    gives the station (metadata) from cadict; match.ghcnst gives the
    station from ghcndict."""

    global wmocount

    for cast in cadict.values():
        result = Struct()
        result.cast = cast
        wmo = cast['WMO Identifier']
        lat = cast['Latitude']
        lon = cast['Longitude']
        caloc = iso6709(lat, lon)
        if wmo:
            result.wmo = wmo
            wmocount += 1
            wmoid11 = "403%s000" % wmo
            if wmoid11 in ghcndict:
                result.ghcnst = ghcndict[wmoid11]
                result.type = 'wmo'
                yield result
                continue
        result.type = 'near'
        nearest = min(ghcndict.values(),
          key=lambda s: locdist((lat,lon), (s.lat, s.lon)))
        result.ghcnst = nearest
        yield result

def drop_short(cadict, cadata):
    """Mutates cadict to drop stations that have fewer than 20 years of
    data."""
    global dropcount
    drops = []
    for key,station in cadict.items():
        id11 = station['id11']
        if id11 not in cadata:
            print "NODATA", id11
            drops.append(key)
            continue
        # Assumes that the scraped data has one duplicate per station.
        data = cadata[id11][id11+'0']
        if len(data) < 20:
            print "SHORT", id11, "YEARS", len(data)
            drops.append(key)
    for station in drops:
        del cadict[station]
    dropcount = len(drops)

def v2asdict(inp):
    """open a v2.mean file and return the data as a dict."""
    def id11(x): return x[:11]
    def id12(x): return x[:12]
    result = {}
    for id11,rows in itertools.groupby(sorted(inp), id11):
        result[id11] = dict((id12, dict((row[12:16],row[16:-1]) for row in rs))
          for id12,rs in itertools.groupby(rows, id12))
    return result

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

def locdist(a, b):
    """*a* and *b* are two locations, each is a (lat,lon) pair in
    decimal degrees.  This function returns the angle between them, in
    decimal degrees (fractional degrees really)."""

    a = map(math.radians, a)
    b = map(math.radians, b)

    def toxyz(p):
        """Convert *p* specified as (lat,lon) in radians to x,y,z
        cartesian."""
        lat,lon = p
        z = math.sin(lat)
        c = math.cos(lat)
        x = c*math.cos(lon)
        y = c*math.sin(lon)
        return x,y,z
    a,b = map(toxyz, [a, b])
    # Use cosine rule to determine angle
    dot = sum(p*q for p,q in zip(a,b))
    # Clamp dot to avoid range error in math.acos due to tiny arithmetic
    # error (dot can be a very little more than 1).
    dot = max(-1, min(1, dot))
    rad = math.acos(dot)
    return math.degrees(rad)

def main(argv=None):
    import sys
    if argv is None:
        argv = sys.argv
    # :todo: replace '403.v2' (which is a grepped selection of v2.mean),
    # with the real v2.mean location
    match(cameta=open('ca.json'), cadata=v2asdict(open('ca.v2')),
      ghcnmeta=open('input/v2.inv'),
      ghcndata=v2asdict(open('403.v2')))

if __name__ == '__main__':
    main()

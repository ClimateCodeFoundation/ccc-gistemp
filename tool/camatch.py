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
import gio

class Struct:
    pass

wmocount = 0

def match(cameta, cadata, ghcnmeta, ghcndata, table):
    """Analyse and match scraped Canada stations with GHCN stations.
    A renaming table is written to *table*.
    """

    cadict = cametaasdict(cameta)
    ghcndict = ghcnmetaasdict(ghcnmeta)

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
        match.sep = sep
        if match.type == 'wmo':
            wmomatch += 1
            if sep > 0.015:
                wmofar += 1
        else:
            if sep < 0.015:
                locnear += 1
        castid = cast['id11']
        overlap,q,id12 = match_quality(cadata[castid][castid+'0'],
            ghcndata[ghcnst.uid])
        match.q = q
        wmo = blah_get(match.cast, 'WMO Identifier') or 'nowmo'
        print match.type, wmo, castid, id12, \
          "%.2f" % sep, "%4d" % overlap, q
        if match.type == 'wmo' or match.q + match.sep < 1:
            newid = ghcnst.uid+'9'
            assert newid not in ghcndata[ghcnst.uid]
            table.write("%s %s\n" % (castid+'0', newid))

    print "dropped stations", dropcount
    print "kept stations", len(cadict)
    print "WMO stations", wmocount
    print "WMO match", wmomatch, "(of those) not near:", wmofar
    print "LOCATION near", locnear

def match_quality(reference, candidate):
    """See how closely the temperature data in *reference* matches the
    temperature data in *candidate*.  *reference* is a single dict
    corresponding to one series; *candidate* is a dict that contains
    multiple ("duplicate") series."""

    ref = as_monthly(reference)
    def scores():
        for id12,record in candidate.items():
            rec = as_monthly(record)
            overlap, q = score(ref, rec)
            yield overlap, q, id12
    def keyscore(t):
        """Takes a tuple *t* (yielded by scores(), above) and converts it
        to a key for comparison.  For long overlaps, we want to score by
        closeness of match (q), for short overlaps, we simply use the
        length of overlap."""

        l,q,id12 = t
        l = min(120, l)
        return l,-q,id12
    return max(scores(), key=keyscore)

def score(ref, cand):
    """Score how close the two series are."""
    overlap = set(ref) & set(cand)
    if not overlap:
        return 0, 9e9
    s = sum((ref[k] - cand[k])**2 for k in overlap)
    q = math.sqrt(s / len(overlap))
    return len(overlap), q

def as_monthly(d):
    """Let *d* be a dict of temperature data such that d[year] is a
    string corresponding to the 12 months of data.  Return a new dict
    *r* where r['YYYY-MM'] gives the value for that month."""

    r = {}
    for year,row in d.items():
        for i in range(12):
            s = row[i*5:(i+1)*5]
            if s == '-9999':
                continue
            x = float(s)*0.1
            m = i+1
            key = "%s-%02d" % (year, m)
            r[key] = x
    return r



def itermatches(cadict, ghcndict):
    """Iterate over all the stations in cadict, yielding a match object
    for each one.  match.type is either 'wmo' or 'loc'; match.cast
    gives the station (metadata) from cadict; match.ghcnst gives the
    station from ghcndict."""

    global wmocount

    for cast in cadict.values():
        result = Struct()
        result.cast = cast
        wmo = blah_get(cast, 'WMO Identifier')
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
        result.type = 'loc'
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
            drops.append((key, 0))
            continue
        # Assumes that the scraped data has one duplicate per station.
        data = cadata[id11][id11+'0']
        if len(data) < 20:
            drops.append((key, len(data)))
    for station,_ in drops:
        del cadict[station]
    print drops
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
            # Convert some key keys into floating point.
            for key in ['Latitude', 'Longitude', 'Elevation']:
                item[key] = float(item[key])
            id = blah_get(item, 'Climate Identifier')
            yield id, item
    return dict(itermeta())


def ghcnmetaasdict(ghcnmeta):
    """Take an open file descriptor on the 'v2.inv' file and return a
    dictionary."""

    stations = gio.station_metadata(file=ghcnmeta, format='v2')
    return stations

def blah_get(target, key):
    """Look up *key* in the dictionary *target*, trying variants of the
    key with spaces replaced with '_'.  Raises exception if a value is
    not found at any of the variants.
    """

    l = [key, key.replace(' ', '_')]
    # Sort any valid keys to the beginning of the list.
    l.sort(key=lambda k: k in target, reverse=True)
    return target[l[0]]

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

def keep403a431(input):
    for row in input:
        if row[:3] in ('403','431'):
            yield row

def main(argv=None):
    import sys
    if argv is None:
        argv = sys.argv
    match(cameta=open('input/ca.json'),
      cadata=v2asdict(open('input/ca.v2.mean')),
      ghcnmeta=open('input/v2.inv'),
      ghcndata=v2asdict(keep403a431(open('input/v2.mean'))),
      table=open('input/ca.tbl', 'w'))

if __name__ == '__main__':
    main()

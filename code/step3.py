#!/usr/bin/env python
# $URL$
# $Rev$
#
# step3.py
#
# David Jones, Ravenbrook Limited, 2008-08-06

"""
Python code reproducing the STEP3 part of the GISTEMP algorithm.
"""

import eqarea
import giss_data
import parameters
import series
from giss_data import MISSING, valid, invalid

import math
# http://docs.python.org/release/2.4.4/lib/module-os.path.html
import os.path
import sys
import itertools

log = open(os.path.join('log', 'step3.log'), 'w')


def incircle(iterable, arc, lat, lon):
    """An iterator that filters iterable (the argument) and yields every
    station with a certain distance of the point of interest given by
    lat and lon (in degrees).  Each station returned has an associated
    weight (normalized distance from centre).  A series of (*station*,
    *weight*) pairs is yielded.

    This is essentially a filter; the stations that are returned are in
    the same order in which they appear in iterable.

    A station record is returned if the great circle arc between it
    and the point of interest is less than *arc* radians (using angles
    makes it independent of sphere size).

    The weight is 1-(d/arc).  where *d* is the
    chord length on a unit circle (from the point lat,lon to the
    station).
    """

    # Warning: lat,lon in degrees; arc in radians!

    cosarc = math.cos(arc)
    coslat = math.cos(lat*math.pi/180)
    sinlat = math.sin(lat*math.pi/180)
    coslon = math.cos(lon*math.pi/180)
    sinlon = math.sin(lon*math.pi/180)

    for record in iterable:
        st = record.station
        s_lat, s_lon = st.lat, st.lon
        # A possible improvement in speed (which the corresponding
        # Fortran code does) would be to store the trig values of
        # the station location in the station object.
        sinlats = math.sin(s_lat*math.pi/180)
        coslats = math.cos(s_lat*math.pi/180)
        sinlons = math.sin(s_lon*math.pi/180)
        coslons = math.cos(s_lon*math.pi/180)

        # Todo: instead of calculating coslon, sinlon, sinlons and coslons,
        # could calculate cos(s_lon - lon),
        # because cosd is (slat1* slat2 + clat1 * clat2*cos(londiff))

        # Cosine of angle subtended by arc between 2 points on a
        # unit sphere is the vector dot product.
        cosd = (sinlats*sinlat +
            coslats*coslat*(coslons*coslon + sinlons*sinlon))
        if cosd > cosarc:
            d = math.sqrt(2*(1-cosd)) # chord length on unit sphere
            weight = 1.0 - (d / arc)
            yield record, weight


def sort(l, cmp):
    """Sort the list l (in place) according to the comparison function
    cmp.  The comparison function, cmp(x, y), should return something
    less than 0 when x < y, 0 when x == y, something greater than 0 when
    x > y.  The sort is ascending in the following sense:  For the sorted
    list: When i < j, cmp(l[i], l[j]) <= 0.

    This sort is not stable.  In fact it is a deliberate emulation of
    the O(n**2) sort implemented in the SORT subroutine of
    to.SBBXgrid.f.  This is necessary to achieve results that are as
    close as possible to the GISS code.  We should switch to using
    Python's built-in sort routine.
    """

    # See to.SBBXgrid.f lines 605 and following

    for n in range(len(l)-1):
        nlmax = n
        for nn in range(n+1, len(l)):
            if cmp(l[nn], l[nlmax]) < 0:
                nlmax = nn
        # swap items at n and nlmax
        t = l[nlmax]
        l[nlmax] = l[n]
        l[n] = t
    return


def iter_subbox_grid(station_records, max_months, first_year, radius):
    """Convert the input *station_records*, into a gridded anomaly
    dataset which is returned as an iterator.

    *max_months* is the maximum number of months in any station
    record.  *first_year* is the first year in the dataset.  *radius*
    is the combining radius in kilometres.
    """

    # Clear Climate Code
    import earth # required for radius.

    # Convert to list because we re-use it for each box (region).
    station_records = list(station_records)
    # Descending sort by number of good records.
    # TODO: Switch to using Python's sort method here, although it
    # will change the results.
    sort(station_records, lambda x,y: y.good_count - x.good_count)

    # A dribble of progress messages.
    dribble = sys.stdout

    # Critical radius as an angle of arc
    arc = radius / earth.radius
    arcdeg = arc * 180 / math.pi

    regions = list(eqarea.gridsub())
    for region in regions:
        box, subboxes = region[0], list(region[1])

        # Count how many cells are empty
        n_empty_cells = 0
        for subbox in subboxes:
            # Select and weight stations
            centre = eqarea.centre(subbox)
            dribble.write("\rsubbox at %+05.1f%+06.1f (%d empty)" % (
              centre + (n_empty_cells,)))
            dribble.flush()
            # Determine the contributing stations to this grid cell.
            contributors = list(incircle(station_records, arc, *centre))

            # Combine data.
            subbox_series = [MISSING] * max_months

            if not contributors:
                box_obj = giss_data.Series(series=subbox_series,
                    box=list(subbox), stations=0, station_months=0,
                    d=MISSING)
                n_empty_cells += 1
                yield box_obj
                continue

            # Initialise series and weight arrays with first station.
            record,wt = contributors[0]
            total_good_months = record.good_count
            total_stations = 1

            offset = record.rel_first_month - 1
            a = record.series # just a temporary
            subbox_series[offset:offset + len(a)] = a
            max_weight = wt
            weight = [wt*valid(v) for v in subbox_series]

            # For logging, keep a list of stations that contributed.
            # Each item in this list is a triple (in list form, so that
            # it can be converted to JSON easily) of [id12, weight,
            # months].  *id12* is the 12 character station identifier;
            # *weight* (a float) is the weight (computed based on
            # distance) of the station's series; *months* is a 12 digit
            # string that records whether each of the 12 months is used.
            # '0' in position *i* indicates that the month was not used,
            # a '1' indicates that is was used.  January is position 0.
            l = [any(valid(v) for v in subbox_series[i::12])
              for i in range(12)]
            s = ''.join('01'[x] for x in l)
            contributed = [[record.uid,wt,s]]

            # Add in the remaining stations
            for record,wt in contributors[1:]:
                # TODO: A method to produce a padded data series
                #       would be good here. Hence we could just do:
                #           new = record.padded_series(max_months)
                new = [MISSING] * max_months
                aa, bb = record.rel_first_month, record.rel_last_month
                new[aa - 1:bb] = record.series
                station_months = series.combine(
                    subbox_series, weight, new, wt,
                    parameters.gridding_min_overlap)
                n_good_months = sum(station_months)
                total_good_months += n_good_months
                if n_good_months == 0:
                    contributed.append([record.uid, 0.0, '0'*12])
                    continue
                total_stations += 1
                s = ''.join('01'[bool(x)] for x in station_months)
                contributed.append([record.uid,wt,s])

                max_weight = max(max_weight, wt)

            series.anomalize(subbox_series,
                             parameters.gridding_reference_period, first_year)
            box_obj = giss_data.Series(series=subbox_series, n=max_months,
                    box=list(subbox), stations=total_stations,
                    station_months=total_good_months,
                    d=radius*(1-max_weight))
            log.write("%s stations %s\n" % (box_obj.uid,
              asjson(contributed)))
            yield box_obj
        plural_suffix = 's'
        if n_empty_cells == 1:
            plural_suffix = ''
        dribble.write(
          '\rRegion (%+03.0f/%+03.0f S/N %+04.0f/%+04.0f W/E): %d empty cell%s.\n' %
            (tuple(box) + (n_empty_cells,plural_suffix)))
    dribble.write("\n")

def asjson(obj):
    """Return a string: The JSON representation of the object "obj".
    This is a peasant's version, not intentended to be fully JSON
    general."""

    return repr(obj).replace("'", '"')


def step3(records, radius=parameters.gridding_radius, year_begin=1880):
    """Step 3 of the GISS processing.

    *records* should be a generator that yields each station.

    """

    # Most of the metadata here used to be synthesized in step2.py and
    # copied from the first yielded record.  Now we synthesize here
    # instead.
    last_year = giss_data.get_last_year()
    year_begin = giss_data.BASE_YEAR
    assert year_begin <= last_year
    # Compute total number of months in a fixed length record.
    monm = 12 * (last_year - year_begin + 1)
    meta = giss_data.SubboxMetaData(mo1=None, kq=1, mavg=6, monm=monm,
            monm4=monm + 7, yrbeg=year_begin, missing_flag=9999,
            precipitation_flag=9999,
            title='GHCN V2 Temperatures (.1 C)')

    units = '(C)'
    title = "%20.20s ANOM %-4s CR %4dKM %s-present" % (meta.title,
            units, radius, year_begin)
    meta.mo1 = 1
    meta.title = title.ljust(80)
    meta.gridding_radius = radius

    box_source = iter_subbox_grid(records, monm, year_begin, radius)

    yield meta
    for box in box_source:
        yield box

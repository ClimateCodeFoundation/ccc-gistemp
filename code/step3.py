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

import earth # required for radius.
import eqarea
import giss_data
import parameters

import math
import sys
import itertools


def invalid(x):
    """Test for invalid datum ("equal" to the giss_data.XMISSING value, for some
    definition of "equal").
    """

    return x == giss_data.XMISSING

def valid(x):
    """Test for valid datum.  See invalid()."""

    # It's important that this obey Iverson's convention: in other words
    # return 0 or 1 (or False or True, which it does).
    return not invalid(x)


def inbox(station_records, lats, latn, longw, longe):
    """An iterator that yields the records for every station within the box
    bounded by the lines of latitude lats (to the south), latn (to the
    north), and the meridians at longw (to the west), and longe (to the
    east).

    In order to accommodate boxes that overlap the meridian at -180
    it is permissible for either longw to be < -180 or for longe to
    be > +180.

    For stations exactly on the boundary the
    "lower-left" rule is used.  Stations are returned if they lie on
    the southern boundary or the western boundary except for
    corners; stations lying exactly on a corner are only returned if
    it is the south-west corner.  Note therefore that to include a
    station situated exactly at the North Pole a latn value slightly
    larger than 90 should be used.
    """

    assert lats <= latn
    assert longw <= longe

    for record in station_records:
        st = record.station
        lat = st.lat
        lon = st.lon

        # if longitude outside box, try mod 360
        if lon > longe:
            lon -= 360
        elif lon < longw:
            lon += 360

        if (lats < lat < latn) and (longw < lon < longe):
            yield record
        if lats == lat and longw <= lon < longe: # southern edge
            yield record
        if longw == lon and lats <= lat < latn: # western edge
            yield record


def incircle(iterable, arc, lat, lon):
    """An iterator that filters iterable (the argument) and yields every
    station with a certain distance of the point of interest given by
    lat and lon (in degrees).

    This is essentially a filter; the stations that are returned are in
    the same order in which they appear in iterable.

    A station record is returned if the great circle arc between it
    and the point of interest is less than arc radians (using angles
    is slightly odd, but makes it independent of sphere size).

    The records returned are given weights (attribute .weight) based
    on its distance from the point.  The weight is 1-(d/arc).  d is 
    not the angle from the station to the point of interest but is the
    chord length on a unit circle.
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
            d /= arc
            record.weight = 1 - d
            yield record


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


def combine(average, weight, new, new_weight, first_year, last_year):
    """Run the GISTEMP combining algorithm.  This combines the data in
    the *new* array into the *average* array.  *new* has weight
    *new_weight*, *average* has weights in the *weight* array.

    Only data for years in *range(first_year, last_year)* are
    considered and combined.

    The number of month records combined is returned.
    
    Each month of the year is considered separately.  For the set of
    times where both *average* and *new* have data the mean difference (a
    bias) is computed.  If there are fewer than
    *parameters.gridding_min_overlap* years in common the data (for
    that month of the year) are not combined.  The bias is subtracted
    from the new record and it is point-wise combined into *average*
    according to the station weight, *new_weight*, and the existing
    weights for average.
    """

    months_combined = 0
    for m in range(12):
        sum_new = 0.0  # Sum of data in new
        sum = 0.0      # Sum of data in average
        count = 0      # Number of years where both new and average are valid
        for a,n in itertools.izip(average[first_year*12+m: last_year*12: 12],
                                  new[first_year*12+m: last_year*12: 12]):
            if invalid(a) or invalid(n):
                continue
            count += 1
            sum += a
            sum_new += n
        if count < parameters.gridding_min_overlap:
            continue
        bias = (sum-sum_new)/count
        # Update period of valid data, averages and weights
        for i in range(first_year*12+m, last_year*12, 12):
            if invalid(new[i]):
                continue
            new_month_weight = weight[i] + new_weight
            average[i] = (weight[i]*average[i] + new_weight*(new[i]+bias))/new_month_weight
            weight[i] = new_month_weight
            months_combined += 1
    return months_combined


# TODO: This was an almost duplicate of the function in step5.
#       Make the code common.
def tavg(data, base, limit):
    """tavg computes the time averages (separately for each calendar
    month) over the base period (year base to limit). In case of no
    data, the average is computed over the whole period.
    """
    averages = []
    for m in range(12):
        sum = 0.0
        count = 0
        for datum in data[m+12*base:m+12*limit:12]:
            if invalid(datum):
                continue
            count += 1
            sum += datum
        if count == 0:
            sum = 0.0
            count = 0
            for datum in data[m::12]:
                if invalid(datum):
                    continue
                count += 1
                sum += datum
        if count == 0:
            averages.append(0.0)
        else:
            averages.append(sum/count)
    return averages


def iter_subbox_grid(station_records, max_months, first_year, radius):
    """Convert the input *station_records*, into a gridded anomaly
    dataset which is returned as an iterator.

    *max_months* is the maximum number of months in any station
     record.  *first_year* is the first year in the dataset.  *radius*
     is the combining radius in kilometres.
    """

    log = sys.stdout

    # Critical radius as an angle of arc
    arc = radius / earth.radius
    arcdeg = arc * 180 / math.pi

    regions = list(eqarea.gridsub())
    for region in regions:
        box, subboxes = region[0], list(region[1])

        # Extend box, by half a box east and west and by arc north
        # and south.
        extent = [box[0] - arcdeg,
                  box[1] + arcdeg,
                  box[2] - 0.5 * (box[3] - box[2]),
                  box[3] + 0.5 * (box[3] - box[2])]
        if box[0] <= -90 or box[1] >= 90:
            # polar
            extent[2] = -180.0
            extent[3] = +180.0

        region_records = list(inbox(station_records, *extent))
        # Descending sort by number of good records
        # TODO: Switch to using Python's sort method here, although it
        # will change the results.
        sort(region_records, lambda x,y: y.good_count - x.good_count)

        # Used to generate the "subbox at" rows in the log.
        lastcentre = (None, None)
        for subbox in subboxes:
            # Convert latitude longitude to integer 100ths for output.
            latlon = map(lambda x:int(round(100*x)), subbox)
            # Select and weight stations
            centre = eqarea.centre(subbox)
            if centre[0] != lastcentre[0]:
                log.write("\nsubbox at %+05.1f" % centre[0])
            log.write('%+06.1f' % centre[1])
            log.flush()
            lastcentre = centre
            # Of possible station records for this region, filter for those
            # from stations within radius of subbox centre.
            incircle_records = list(incircle(region_records, arc, *centre))

            # Combine data.
            series = [giss_data.XMISSING] * max_months

            if len(incircle_records) == 0:
                box_obj = giss_data.SubboxRecord(
                    lat_S=latlon[0], lat_N=latlon[1], lon_W=latlon[2],
                    lon_E=latlon[3], stations=0, station_months=0,
                    d=giss_data.XMISSING, series=series)
                log.write('*')
                log.flush()
                yield box_obj
                continue

            # Initialise data with first station
            record = incircle_records[0]
            total_good_months = record.good_count
            total_stations = 1

            max_weight = record.weight
            offset = record.rel_first_month - 1
            a = record.series # just a temporary
            series[offset:offset + len(a)] = a
            weight = [0.0] * max_months
            for i in range(len(a)):
                if valid(a[i]):
                    weight[i + offset] = record.weight

            # Add in the remaining stations
            for record in incircle_records[1:]:
                # TODO: A StationMethod method to produce a padded data series
                #       would be good here. Hence we could just do:
                #           new = record.padded_series(max_months)
                new = [giss_data.XMISSING] * max_months
                aa, bb = record.rel_first_month, record.rel_last_month
                new[aa - 1:bb] = record.series
                station_months = combine(series, weight, new, record.weight, record.rel_first_year, record.rel_last_year + 1)
                total_good_months += station_months
                if station_months == 0:
                    continue
                total_stations += 1

                if max_weight < record.weight:
                    max_weight = record.weight

            bias = tavg(series, parameters.gridding_reference_first_year-first_year,
                        parameters.gridding_reference_last_year-first_year+1)
            m = 0
            for y in range(max_months // 12):
                for k in range(12):
                    if valid(series[m]):
                        series[m] -= bias[k]
                    m += 1
            box_obj = giss_data.SubboxRecord(n=max_months,
                    lat_S=latlon[0], lat_N=latlon[1], lon_W=latlon[2],
                    lon_E=latlon[3], stations=total_stations, station_months=total_good_months,
                    d=radius*(1-max_weight), series=series)
            yield box_obj

    print >>log


def step3(record_source, radius=parameters.gridding_radius, year_begin=1880):
    """Step 3 of the GISS processing.

    """
    station_records = [record for record in record_source]
    m, station_records = station_records[0], station_records[1:]
    meta = giss_data.SubboxMetaData(m.mo1, m.kq, m.mavg, m.monm,
            m.monm + 7, m.yrbeg, m.missing_flag, m.precipitation_flag,
            m.title)

    # Output series cannot begin earlier than input series.
    # TODO: The ``year_begin`` argument seems to have no effect,
    #       other than what appears in the ``title`` below.
    #       This might have been inroduced during Jan/Feb 2010.
    year_begin = max(year_begin, meta.yrbeg)

    units = '(C)'
    title = "%20.20s ANOM %-4s CR %4dKM %s-present" % (meta.title,
            units, radius, year_begin)
    meta.mo1 = 1
    meta.title = title.ljust(80)

    box_source = iter_subbox_grid(station_records,
                                  meta.monm, meta.yrbeg,
                                  radius)

    yield meta
    for box in box_source:
        yield box

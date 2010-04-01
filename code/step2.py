# $URL$
# $Rev$
#
# step2.py
#
# Clear Climate Code, 2010-01-22

"""Python replacement for code/STEP2/*
"""
__docformat__ = "restructuredtext"


import sys
import math
import itertools

import earth
import giss_data
import parameters
from giss_data import valid, invalid, MISSING

def drop_short_records(record_source):
    """Drop records which don't have at least one month index
    with at least *parameters.station_drop_minimum_months* valid data.
    """
    for record in record_source:
        mmax = max(record.get_monthly_valid_counts())
        if mmax >= parameters.station_drop_minimum_months:
            yield record

def annual_anomaly(record):
    """Updates the station record *record* with its annual anomalies.  The
    .anomalies attribute is assigned a list of annual anomalies.
    The .first attribute is assigned the year to which to first item in
    the .anomalies list applies.

    The algorithm is as follows: compute monthly averages, then
    monthly anomalies, then seasonal anomalies (means of monthly
    anomalies for at least two months) then annual anomalies (means
    of seasonal anomalies for at least three seasons).
    """

    series = record.series
    monthly_means = []
    for m in range(12):
        month_data = filter(valid, series[m::12])
        # neglect December of final year, as we do not use its season.
        if m == 11 and valid(series[-1]):
            month_data = month_data[:-1]
        monthly_means.append(float(sum(month_data)) / len(month_data))
    annual_anoms = []
    first = None
    for y in range(len(series)/12):
        # Seasons are Dec-Feb, Mar-May, Jun-Aug, Sep-Nov.
        # (Dec from previous year).
        total = [0.0] * 4 # total monthly anomaly for each season
        count = [0] * 4   # number of valid months in each season
        for m in range(-1, 11):
            index = y * 12 + m
            if index >= 0: # no Dec value in year -1
                datum = series[index]
                if valid(datum):
                    season = (m+1) // 3 # season number 0-3
                    total[season] += datum - monthly_means[(m + 12) % 12]
                    count[season] += 1
        season_anomalies = [] # list of valid seasonal anomalies
        for s in range(4):
            # valid seasonal anomaly requires at least 2 valid months
            if count[s] > 1:
                season_anomalies.append(total[s]/count[s])
        # valid annual anomaly requires at least 3 valid seasons
        if len(season_anomalies) > 2:
            annual_anoms.append(sum(season_anomalies) / len(season_anomalies))
            if first is None:
                first = y
            last = y
        else:
            annual_anoms.append(MISSING)

    if first is None:
        record.anomalies = None
    else:
        record.first = first + record.first_year
        record.anomalies = annual_anoms[first: last+1]


def is_rural(station):
    """Test whether the station described by *station* is rural.
    """
    if parameters.use_global_brightness:
        return station.global_brightness <= 10
    else:
        return (station.US_brightness == '1' or
                (station.US_brightness == ' ' and
                 station.pop == 'R'))

class Struct(object):
    pass

def urban_adjustments(anomaly_stream):
    """Takes an iterator of station records and applies an adjustment
    to urban stations to compensate for urban temperature effects.
    Returns an iterator of station records.  Rural stations are passed
    unchanged.  Urban stations which cannot be adjusted are discarded.

    The adjustment follows a linear or two-part linear fit to the
    difference in annual anomalies between the urban station and the
    combined set of nearby rural stations.  The linear fit is to allow
    for a linear effect at the urban station.  The two-part linear fit
    is to allow for a model of urban effect which starts or stops at
    some point during the time series.

    The algorithm is essentially as follows:

    For each urban station:
        1. Find all the rural stations within a fixed radius;
        2. Combine the annual anomaly series for those rural stations, in
           order of valid-data count;
        3. Calculate a two-part linear fit for the difference between
           the urban annual anomalies and this combined rural annual anomaly;
        4. If this fit is satisfactory, apply it; otherwise apply a linear fit.

        If there are not enough nearby rural stations, or the combined
        rural record does not have enough overlap with the urban
        record, try a second time for this urban station, with a
        larger radius.  If there is still not enough data, discard the
        urban station.
     """

    last_year = giss_data.get_ghcn_last_year()
    first_year = 1880

    iyoff = giss_data.BASE_YEAR - 1
    iyrm = last_year - iyoff

    rural_stations = []
    urban_stations = {}

    pi180 = math.pi / 180.0

    all = []
    for record in anomaly_stream:
        station = record.station
        all.append(record)
        record.urban_adjustment = None
        annual_anomaly(record)
        if record.anomalies is None:
            continue
        length = len(record.anomalies)
        d = Struct()
        d.anomalies = record.anomalies
        d.cslat = math.cos(station.lat * pi180)
        d.snlat = math.sin(station.lat * pi180)
        d.cslon = math.cos(station.lon * pi180)
        d.snlon = math.sin(station.lon * pi180)
        d.id = record.uid
        d.first_year = record.first - iyoff
        d.last_year = d.first_year + length - 1
        d.station = station
        d.record = record
        if is_rural(station):
            rural_stations.append(d)
        else:
            urban_stations[record] = d

    # Sort the rural stations according to the length of the time record
    # (ignoring gaps).
    for st in rural_stations:
        st.recLen = len([v for v in st.anomalies if valid(v)])
    rural_stations.sort(key=lambda s:s.recLen)
    rural_stations.reverse()

    # Combine time series for rural stations around each urban station
    for record in all:
        us = urban_stations.get(record, None)
        if us is None:
            # Just remove leading/trailing invalid values for rural stations.
            record.strip_invalid()
            yield record
            continue

        iyu1 = us.first_year + iyoff - 1 # subtract 1 for a possible partial yr
        iyu2 = us.last_year + iyoff + 1  # add 1 for partial year

        usingFullRadius = False
        dropStation = False
        needNewNeighbours = True
        while True:
            if needNewNeighbours:
                if usingFullRadius:
                    radius = parameters.urban_adjustment_full_radius
                else:
                    radius = parameters.urban_adjustment_full_radius / 2
                neighbors = get_neighbours(us, rural_stations, radius)
                if not neighbors:
                    if usingFullRadius:
                        dropStation = True
                        break
                    usingFullRadius = True
                    needNewNeighbours = True
                    continue

                counts, urban_series, combined = combine_neighbors(
                        us, iyrm, iyoff, neighbors)
                iy1 = 1
                needNewNeighbours = False

            points, quorate_count, first, last = prepare_series(
                iy1, iyrm, combined, urban_series, counts, iyoff)

            if quorate_count < parameters.urban_adjustment_min_years:
                if usingFullRadius:
                    dropStation = True
                    break
                usingFullRadius = True
                needNewNeighbours = True
                continue

            if quorate_count >= (parameters.urban_adjustment_proportion_good
                                 * (last - first + 0.9)):
                break

            # Not enough good years for the given range.  Try to save
            # cases in which the gaps are in the early part, by
            # dropping that part and going around to prepare_series
            # again.
            iy1 = int(last - (quorate_count - 1) /
                      parameters.urban_adjustment_proportion_good)
            if iy1 < first + 1:
                iy1 = first + 1                  # avoid infinite loop

        if dropStation:
            continue

        fit = getfit(points)
        # find extended range
        iyxtnd = int(round(quorate_count /
                           parameters.urban_adjustment_proportion_good)
                     - (last - first + 1))
        n1x = first + iyoff
        n2x = last + iyoff
        if iyxtnd < 0:
            sys.exit('impossible')
        if iyxtnd > 0:
            lxend = iyu2 - (last + iyoff)
            if iyxtnd <= lxend:
                 n2x = n2x + lxend
            else:
                 n1x = n1x - (iyxtnd - lxend)
                 if n1x < iyu1:
                     n1x = iyu1
                 n2x = iyu2

        series = record.series
        # adjust
        m1 = record.rel_first_month + record.good_start_idx
        m2 = record.rel_first_month + record.good_end_idx - 1
        offset = record.good_start_idx # index of first valid month
        a, b = adjust(first_year, record, series, fit, n1x, n2x,
                      first + iyoff, last + iyoff, m1, m2, offset)
        # a and b are numbers of new first and last valid months
        aa = a - m1
        bb = b - a + 1
        record.set_series(a-1 + first_year * 12 + 1,
                          series[aa + offset:aa + offset + bb])
        del record.first
        yield record


def get_neighbours(us, rural_stations, radius):
    """Returns a list of the stations in *rural_stations* which are
    within distance *radius* of the urban station *us*.  Each rural
    station returned is given a 'weight' slot representing its
    distance fromn the urban station.
    """
    neighbors = []

    cos_crit = math.cos(radius / earth.radius)
    rbyrc = earth.radius / radius

    for rs in rural_stations:
        csdbyr = (rs.snlat * us.snlat + rs.cslat * us.cslat *
                     (rs.cslon * us.cslon  + rs.snlon * us.snlon))

        if csdbyr <= cos_crit:
            continue
        dbyrc = 0
        if csdbyr < 1.0:
            dbyrc = rbyrc * math.sqrt(2.0 * (1.0 - csdbyr))
        rs.weight = 1.0 - dbyrc
        neighbors.append(rs)

    return neighbors


def combine_neighbors(us, iyrm, iyoff, neighbors):
    """Combines the neighbor stations *neighbors*, weighted according
    to their distances from the urban station *us*, to give a combined
    annual anomaly series.  Returns a tuple: (*counts*,
    *urban_series*, *combined*), where *counts* is a per-year list of
    the number of stations combined, *urban_series* is the series from
    the urban station, re-based at *iyoff*, and *combined* is the
    combined neighbor series, based at *iyoff*.
    """

    weights = [0.0] * iyrm
    counts = [0] * iyrm
    urban_series = [MISSING] * iyrm
    combined = [MISSING] * iyrm

    urban_series[us.first_year - 1:us.last_year] = us.anomalies

    # start with the neighbor with the longest time record
    rs = neighbors[0]
    combined[rs.first_year - 1:rs.last_year] = rs.anomalies
    for m in range(len(rs.anomalies)):
        if valid(rs.anomalies[m]):
            weights[m + rs.first_year - 1] = rs.weight
            counts[m + rs.first_year - 1] = 1

    # add in the remaining stations
    for i, rs in enumerate(neighbors[1:]):
        dnew = [MISSING] * iyrm
        dnew[rs.first_year - 1: rs.last_year] = rs.anomalies
        cmbine(combined, weights, counts, dnew,
               rs.first_year, rs.last_year, rs.weight)

    return counts, urban_series, combined


def prepare_series(iy1, iyrm, combined, urban_series, counts, iyoff):
    """Prepares for the linearity fitting by returning a series of
    data points *(x,f)*, where *x* is a year number and *f* is the
    difference between the combined rural station anomaly series
    *combined* and the urban station series *urban_series*.  The
    points only include valid years, from the first quorate year to
    the last.  A valid year is one in which both the urban station and
    the combined rural series have valid data.  A quorate year is a
    valid year in which there are at least
    *parameters.urban_adjustment_min_rural_stations* contributing
    (obtained from the *counts* series).

    Returns a 4-tuple: (*p*, *c*, *f*, *l*). *p* is the series of
    points, *c* is a count of the valid quorate years.  *f* is the
    first such year.  *l* is the last such year.
    """
    first = last = i = quorate_count = length = 0
    points = []

    for iy in xrange(iy1 - 1, iyrm):
        if valid(combined[iy]) and valid(urban_series[iy]):
            if counts[iy] >= parameters.urban_adjustment_min_rural_stations:
                last = iy + 1
                quorate_count += 1
                if first == 0:
                    first = iy + 1
            if quorate_count <= 0:
                continue

            points.append((iy + iyoff + 1, combined[iy] - urban_series[iy]))
            i += 1
            if counts[iy] >= parameters.urban_adjustment_min_rural_stations:
                 length = i

    return points[:length], quorate_count, first, last


def cmbine(combined, weights, counts, data, first, last, weight):
    """Adds the array *data* with weight *weight* into the array of
    weighted averages *combined*, with total weights *weights* and
    combined counts *counts* (that is, entry *combined[i]* is the
    result of combining *counts[i]* values with total weights
    *weights[i]*).  Adds the computed bias between *combined* and
    *data* before combining.

    Only combines in the range [*first*, *last*); only combines valid
    values from *data*, and if there are fewer than
    *parameters.rural_station_min_overlap* entries valid in both
    arrays then it doesn't combine at all.

    Note: if *data[i]* is valid and *combined[i]* is not, the weighted
    average code runs and still produces the right answer, because
    *weights[i]* will be zero.
    """
    sumn = ncom = 0
    avg_sum = 0.0
    a, b = first - 1, last
    for v_avg, v_new in itertools.izip(combined[a:b], data[a:b]):
        if invalid(v_avg) or invalid(v_new):
            continue
        ncom = ncom + 1
        avg_sum += v_avg
        sumn += v_new

    if ncom < parameters.rural_station_min_overlap:
        return
    bias = (avg_sum - sumn) / float(ncom)

    # update period of valid data, averages and weights
    for n in xrange(first - 1, last):
        v_new = data[n]
        if invalid(v_new):
            continue
        wtnew = weights[n] + weight
        old_wt, weights[n] = weights[n], wtnew
        combined[n] = (old_wt * combined[n] + weight * (v_new + bias)) / wtnew
        counts[n] += 1


def getfit(points):
    """ Finds the best two-part linear fit for *points*, and
    returns the fit parameters.
    """

    # Todo: incorporate trend2 into this.
    rmsmin = 1.e20

    for n in xrange(parameters.urban_adjustment_min_leg,
                    len(points) - parameters.urban_adjustment_min_leg):
        xknee = points[n][0]
        sl1, sl2, rms, sl = trend2(points, xknee, 2)

        if rms < rmsmin:
             rmsmin = rms
             fit = (sl1, sl2, xknee, sl)

    return fit


def trend2(points, xmid, min):
    """Finds a fit to the data *points[]*, using regression analysis,
    by a line with a change in slope at *xmid*. Returned is a 4-tuple
    (*sl1*, *sl2*, *rms*, *sl*): the left-hand slope, the right-hand
    slope, the RMS error, and the slope of an overall linear fit.
    """

    # Todo: incorporate into getfit.
    count0 = count1 = 0
    sx0 = sx1 = 0
    sxx0 = sxx1 = 0
    sxa0 = sxa1 = 0

    sa = 0.0
    saa = 0.0

    for (x,v) in points:
        if invalid(v):
            continue
        x -= xmid
        sa += v
        saa += v ** 2
        if x > 0.0:
            count1 += 1
            sx1 += x
            sxx1 += x ** 2
            sxa1 += x * v
        else:
            count0 += 1
            sx0 += x
            sxx0 += x ** 2
            sxa0 += x * v

    if count0 < min or count1 < min:
       return MISSING, MISSING, MISSING, MISSING

    count = count0 + count1
    denom = (count * sxx0 * sxx1
             - sxx0 * sx1 ** 2
             - sxx1 * sx0 ** 2)
    sl1 = (sx0 * (sx1 * sxa1 - sxx1 * sa)
           + sxa0 * (count * sxx1 - sx1 ** 2)) / denom
    sl2 = (sx1 * (sx0 * sxa0 - sxx0 * sa)
           + sxa1 * (count * sxx0 - sx0 ** 2)) / denom

    ymid = (sa - sl1 * sx0 - sl2 * sx1) / count
    rms = (count * ymid ** 2
           + saa
           - 2 * ymid * (sa - sl1 * sx0 - sl2 * sx1)
           + sl1 * sl1 * sxx0
           + sl2 * sl2 * sxx1
           - 2 * sl1 * sxa0
           - 2 * sl2 * sxa1)

    # linear regression
    sx = sx0 + sx1
    sxx = sxx0 + sxx1
    sxa = sxa0 + sxa1
    sl = (count * sxa - sa * sx) / (count * sxx - sx ** 2)

    return sl1, sl2, rms, sl


def adjust(first_year, station, series, fit, iy1, iy2,
           iy1a, iy2a, m1, m2, offset):
    (sl1, sl2, knee, sl0) = fit
    if not good_two_part_fit(fit, iy1a, iy2a):
        # Use linear approximation
        sl1, sl2 = sl0, sl0

    base = m1

    m1o, m2o = m1, m2
    m1 = -100
    m0 = 12 * (iy1 - first_year)   # Dec of year iy1
    for iy in range(iy1, iy2 + 1):
        sl = sl1
        if iy > knee:
            sl = sl2
        iya = iy
        if iy < iy1a:
            iya = iy1a
        if iy > iy2a:
            iya = iy2a
        adj = (iya - knee) * sl - (iy2a - knee) * sl2
        for m in range(m0, m0 + 12):
            mIdx = m - base
            if mIdx < 0:
                continue
            if m >= m1o and m <= m2o and valid(series[mIdx + offset]):
                if m1 < 0:
                    m1 = m
                series[mIdx+offset] = series[mIdx+offset] + adj
                m2 = m

        m0 = m0 + 12

    return m1, m2


def good_two_part_fit(fit, iy1, iy2):
    """Decide whether to apply a two-part fit.

    If the two-part fit is not good, the linear fit is used instead.
    The two-part fit is good if all of these conditions are true:

    - left leg is longer than urban_adjustment_short_leg
    - right leg is longer than urban_adjustment_short_leg
    - left gradient is abs less than urban_adjustment_steep_leg
    - right gradient is abs less than urban_adjustment_steep_leg
    - difference between gradients is abs less than urban_adjustment_steep_leg
    - either gradients have same sign or
             at least one gradient is abs less than
             urban_adjustment_reverse_gradient
    """

    (sl1, sl2, knee, sl) = fit
    return ((knee >= iy1 + parameters.urban_adjustment_short_leg) and
            (knee <= iy2 - parameters.urban_adjustment_short_leg) and
            (abs(sl1) <= parameters.urban_adjustment_steep_leg) and
            (abs(sl2) <= parameters.urban_adjustment_steep_leg) and
            (abs(sl2 - sl1) <= parameters.urban_adjustment_steep_leg) and
            ((sl1 * sl2 >= 0) or
             (abs(sl1) <= parameters.urban_adjustment_reverse_gradient) or
             (abs(sl2) <= parameters.urban_adjustment_reverse_gradient)))


def step2(record_source):
    data = drop_short_records(record_source)
    adjusted = urban_adjustments(data)
    for record in adjusted:
        yield record

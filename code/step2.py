# $URL$
# $Rev$
#
# step2.py
#
# Clear Climate Code, 2010-01-22

"""Python replacement for code/STEP2/*

This step performs a trend adjustment to those stations identified as
urban.  An urban station must meet various criteria to have its trend
adjusted (there must be sufficient nearby rural stations and their
combined record must have enough overlap with the urban station); if an
urban station does not meet the criteria, it is discarded.

Additionally in this step, all stations with short records are discarded
(prior to any other processing).
"""
__docformat__ = "restructuredtext"


# Standard Python
import math
import itertools

# Clear Climate Code
import earth
import giss_data
import parameters
from giss_data import valid, invalid, MISSING


def urban_adjustments(record_stream):
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

    iyoff = giss_data.BASE_YEAR - 1
    iyrm = last_year - iyoff

    rural_stations = []
    urban_stations = {}

    pi180 = math.pi / 180.0

    all = []
    for record in record_stream:
        station = record.station
        all.append(record)
        anomalies, first = annual_anomaly(record)
        if anomalies is None:
            continue
        d = Struct()
        d.anomalies = anomalies
        d.cslat = math.cos(station.lat * pi180)
        d.snlat = math.sin(station.lat * pi180)
        d.cslon = math.cos(station.lon * pi180)
        d.snlon = math.sin(station.lon * pi180)
        d.id = record.uid
        length = len(d.anomalies)
        # first_year and last_year are the time bounds of the anomaly
        # series, where 1 is 1880 (BASE_YEAR really).
        d.first_year = first - iyoff
        d.last_year = d.first_year + length - 1
        d.record = record
        if is_rural(station):
            rural_stations.append(d)
        else:
            urban_stations[record] = d

    # Sort the rural stations according to the length of the time record
    # (ignoring gaps).
    for st in rural_stations:
        st.recLen = len([v for v in st.anomalies if valid(v)])
    # Note: Changing the following to use `reverse=True` will change the
    # results (a little bit), because the list will end up in a
    # different order.
    rural_stations.sort(key=lambda s:s.recLen)
    rural_stations.reverse()

    # Combine time series for rural stations around each urban station
    for record in all:
        us = urban_stations.get(record, None)
        if us is None:
            # Not an urban station.  Pass through unchanged.
            yield record
            continue

        dropStation = True
        R = parameters.urban_adjustment_full_radius
        for radius in [R/2, R]:
            neighbours = get_neighbours(us, rural_stations, radius)
            if not neighbours:
                continue
            counts, urban_series, combined = combine_neighbours(
                    us, iyrm, neighbours)
            start_year = 1

            while True:
                points, quorate_count, first, last = prepare_series(
                    start_year, combined, urban_series, counts)

                if quorate_count < parameters.urban_adjustment_min_years:
                    break

                if quorate_count >= (parameters.urban_adjustment_proportion_good
                                     * (last - first + 0.9)):
                    break

                # Not enough good years for the given range.  Try to save
                # cases in which the gaps are in the early part, by
                # dropping that part and going around to prepare_series
                # again.
                start_year = int(last - (quorate_count - 1) /
                          parameters.urban_adjustment_proportion_good)
                start_year = max(start_year, first + 1)

            # Now work out why we exited previous loop.
            if quorate_count < parameters.urban_adjustment_min_years:
                continue
            # We have a good adjustment.
            assert quorate_count >= (parameters.urban_adjustment_proportion_good
                                 * (last - first + 0.9))
            dropStation = False
            break

        if dropStation:
            continue

        fit = getfit(points)
        # find extended range
        iyxtnd = int(round(quorate_count /
                           parameters.urban_adjustment_proportion_good)
                     - (last - first + 1))
        assert iyxtnd >= 0
        # *n1x* and *n2x* are calendar years.
        n1x = first + iyoff
        n2x = last + iyoff

        iyu1 = us.first_year + iyoff - 1 # subtract 1 for a possible partial yr
        iyu2 = us.last_year + iyoff + 1  # add 1 for partial year

        if iyxtnd > 0:
            # When extending, extend to include all of the recent part
            # of the urban record...
            lxend = iyu2 - n2x
            if iyxtnd > lxend:
                # ... and if we have enough "spare years" extend some or
                # all of the earlier part of the urban record.
                n1x = n1x - (iyxtnd - lxend)
                n1x = max(n1x, iyu1)
            n2x = iyu2

        series = record.series
        # adjust
        m1 = record.rel_first_month + record.good_start_idx
        m2 = record.rel_first_month + record.good_end_idx - 1
        offset = record.good_start_idx # index of first valid month
        a, b = adjust(series, fit, n1x, n2x,
                      first + iyoff, last + iyoff, m1, m2, offset)
        # *a* and *b* are numbers of new first and last valid months

        # *lpad* is the number of months to remove starting with the
        # first valid month of the series.
        lpad = a - m1
        # *nretain* is the number of months (valid or invalid) to
        # retain.
        nretain = b - a + 1
        # Truncate the series.
        record.set_series(a-1 + giss_data.BASE_YEAR * 12 + 1,
                          series[lpad + offset:lpad + offset + nretain])
        yield record


def annual_anomaly(record):
    """Computes annual anomalies for the station record *record*.
    Returns a pair (*anomalies*, *first*) where *anomalies* is a
    list of annual anomalies, and *first* is the year to which to
    first item in the list applies.

    If no anomalies can be computed, then (None, None) is returned.

    The algorithm is as follows: compute monthly averages, then
    monthly anomalies, then seasonal anomalies (means of monthly
    anomalies for at least two months) then annual anomalies (means
    of seasonal anomalies for at least three seasons).

    This function assumes that the series starts in January.
    """
    # :todo: Would this be simpler if it didn't truncate the returned
    # series, returning a full length series with MISSING at beginning
    # and end?

    series = record.series
    monthly_means = []
    for m in range(12):
        month_data = series[m::12]
        # Neglect December of final year, as we do not use its season.
        if m == 11:
            month_data = month_data[:-1]
        month_data = filter(valid, month_data)
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
                    # season number 0-3
                    season = (m+1) // 3
                    total[season] += datum - monthly_means[m % 12]
                    count[season] += 1
        season_anomalies = [] # list of valid seasonal anomalies
        for s in range(4):
            # valid seasonal anomaly requires at least 2 valid months
            if count[s] >= 2:
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
        return None, None
    else:
        return (annual_anoms[first: last+1], first + record.first_year)


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


def get_neighbours(us, rural_stations, radius):
    """Returns a list of the stations in *rural_stations* which are
    within distance *radius* of the urban station *us*.  Each rural
    station returned is given a 'weight' slot representing its
    distance fromn the urban station.
    """
    neighbours = []

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
        neighbours.append(rs)

    return neighbours


def combine_neighbours(us, iyrm, neighbours):
    """Combines the neighbour stations *neighbours*, weighted according
    to their distances from the urban station *us*, to give a combined
    annual anomaly series.  *iyrm* is the length of the combined series.
    
    Returns a tuple: (*counts*, *urban_series*, *combined*), where
    *counts* is a per-year list of the number of stations combined,
    *urban_series* is the series from the urban station, re-based.
    *combined* is the combined neighbour series, re-based.
    """

    weights = [0.0] * iyrm
    counts = [0] * iyrm
    urban_series = [MISSING] * iyrm
    combined = [MISSING] * iyrm

    urban_series[us.first_year - 1: us.last_year] = us.anomalies

    # start with the neighbour with the longest time record
    rs = neighbours[0]
    combined[rs.first_year - 1:rs.last_year] = rs.anomalies
    for m,anom in enumerate(rs.anomalies):
        if valid(anom):
            weights[m + rs.first_year - 1] = rs.weight
            counts[m + rs.first_year - 1] = 1

    # add in the remaining stations
    for rs in neighbours[1:]:
        dnew = [MISSING] * iyrm
        dnew[rs.first_year - 1: rs.last_year] = rs.anomalies
        cmbine(combined, weights, counts, dnew,
               rs.first_year, rs.last_year, rs.weight)

    return counts, urban_series, combined


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


def prepare_series(iy1, combined, urban_series, counts):
    """Prepares for the linearity fitting by returning a series of
    data points *(x,d)*, where *x* is a calendar year number and *d*
    is the difference between the combined rural station anomaly series
    *combined* and the urban station series *urban_series* (each of
    these is an annual series, one datum per year).
    
    The returned points only include valid years, from the first
    quorate year to the last quorate year.  A valid year is one in
    which both the urban station and the combined rural series have
    valid data.  A quorate year is a valid year in which there are
    at least *parameters.urban_adjustment_min_rural_stations*
    contributing.

    The algorithm is restricted to only considering years starting at
    *iy1* (and ending at the end of the series); *iy1* is a
    1-based index (for historical reasons).

    The *counts* argument is a sequence that contains the number of
    stations contributing to each datum in *combined*.

    Returns a 4-tuple: (*p*, *c*, *f*, *l*). *p* is the series of
    points, *c* is a count of the valid quorate years.  *f* is the
    first such year.  *l* is the last such year.

    For historical reasons *f* and *l* are returned as 1-based indexes.
    """

    # Calendar year corresponding to first datum in series.
    year_offset = giss_data.BASE_YEAR
    first = last = quorate_count = length = 0
    points = []

    assert len(combined) == len(urban_series)

    def quorate():
        """True when *iy* corresponds to a quorate year; used in inner
        loop, below."""
        return counts[iy] >= parameters.urban_adjustment_min_rural_stations

    for iy in xrange(iy1 - 1, len(combined)):
        if valid(combined[iy]) and valid(urban_series[iy]):
            if quorate():
                last = iy + 1
                quorate_count += 1
                if first == 0:
                    first = iy + 1
            if quorate_count == 0:
                continue

            points.append((iy + year_offset, combined[iy] - urban_series[iy]))
            if quorate():
                 length = len(points)

    return points[:length], quorate_count, first, last


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


def adjust(series, fit, iy1, iy2,
           iy1a, iy2a, m1, m2, offset):
    """Adjust the series according to the previously computed
    parameters.
    
    *series* is a monthly data series;  it is mutated, but its
    length is not changed.
    
    *iy1*, *iy2* are calendar years: the first and last years that are
    subject to adjustment.

    *iy1a*, *iy2a* are calendar years: the first and last years for the
    2-part fit.

    *m1*, *m2* are month numbers for the first and last month with good
    data (where 0 is 1880-01).

    *offset* is the index of the first valid month.

    """
    first_year = giss_data.BASE_YEAR

    # *sl1* and *sl2* are the slopes, in Kelvin per year, of the first
    # and second parts of the fit.  *knee* is the calendar year at which
    # the two slopes meet.
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
                series[mIdx+offset] += adj
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


def drop_short_records(record_source):
    """Drop records which don't have at least one month index
    with at least *parameters.station_drop_minimum_months* valid data.
    """
    for record in record_source:
        mmax = max(record.get_monthly_valid_counts())
        if mmax >= parameters.station_drop_minimum_months:
            yield record

def step2(record_source):
    data = drop_short_records(record_source)
    adjusted = urban_adjustments(data)
    for record in adjusted:
        yield record

#!/usr/bin/env python
# $URL: https://ccc-gistemp.googlecode.com/svn/trunk/code/invnt.py $
# $Rev: 162 $
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
import read_config

import earth
import fort
import ccc_binary

BAD = 9999

def invalid(x):
    """Test for invalid datum ("equal" to the BAD value, for some
    definition of "equal").
    """

    return abs(x - BAD) < 0.1

def valid(x):
    """Test for valid datum.  See invalid()."""

    # It's important that this obey Iverson's convention: in other words
    # return 0 or 1 (or False or True, which it does).
    return not invalid(x)

def read_text():
    """ Reads the data file work/Ts.txt, output by step 1, and returns
    an iterator of the contents.

    Each item in the iterator is a pair (dict, series).  dict contains
    station information, and series contains the temperature series,
    as a list of monthly values.  Each value is an integer in tenths
    of degrees C.

    Logs to log/short.station.list and log/station.log as we go; these
    files may not be necessary, and their formats are a little
    awkward, but retained for now for testing against Fortran.
    """

    v2_inv = read_config.v2_get_info()
    in_file = open("work/Ts.txt")
    short_station_log = open('log/short.station.list', "w")
    station_log = open('log/station.log', "w")
    mult = 0
    last_id11 = None
    minimum_monthly_max = 20
    min_year = 1880

    for (station_line, lines) in itertools.groupby(in_file, lambda line: line[0] == ' '):
        if station_line: # line beginning with a blank introduces a new station
            lines = list(lines)
            assert len(lines) == 1
            line = lines[0]
            id12 = line[10:22]
            id11 = id12[:11]
            if id11 == last_id11:
                mult += 1
            else:
                mult = 0
                last_id11 = id11
            dict = v2_inv[id11].copy()
            dict['id'] = id12
        else: # lines consists of the temperature series
            series = []
            monthly_valid = [0] * 12
            min_month = None
            for line in lines:
                data = map(int, line.split())
                year = int(data[0])
                if not dict.has_key('begin'):
                    dict['begin'] = year
                monthlies = data[1:]
                for month in range(12):
                    if valid(monthlies[month]):
                        monthly_valid[month] += 1
                        if min_month is None:
                            min_month = month + 12*year
                        max_month = month + 12*year
                series.extend(monthlies)
            mmax = max(monthly_valid)
            station_log.write("%12d%12d%12d%12d%12d%12d\n" % (
                    mult, int(dict['id'][3:]), mmax, minimum_monthly_max,
                    min_month - min_year * 12 + 1,
                    max_month - min_year * 12 + 1))
            if mmax >= minimum_monthly_max:
                dict['end'] = year
                dict['years'] = dict['end'] - dict['begin'] + 1
                dict['min_month'] = min_month
                dict['max_month'] = max_month
                yield (dict, series)
            else:
                short_station_log.write("%12d  %30s%c%c%c%3s dropped\n" %
                                        (int(dict['id'][3:]), dict['name'],
                                         dict['US-brightness'],
                                         dict['pop'], dict['GHCN-brightness'],
                                         dict['id'][:3]))

def invnt(name, stream):
    """Turns the data *stream* into a station list, in log/*name*.station.list,
    while passing the stream through unchanged.  Probably unnecessary, as the
    .station.list files are not used by any other code, but retained for now for
    testing against Fortran.
    """
    
    log = open('log/%s.station.list' % name,'w')
    for (dict, series) in stream:
        log.write("work/%s %9d %-30s lat,lon (.1deg)%5d%6d %s%s%s cc=%3s\n" % (
                  name, int(dict['id'][3:12]), dict['name'],
                  math.floor(dict['lat'] * 10.0 + 0.5),
                  math.floor(dict['lon'] * 10.0 + 0.5),
                  dict['pop'], dict['US-brightness'], dict['GHCN-brightness'], dict['id'][:3]))
        yield (dict, series)
    log.close()

def annual_anomalies(stream):
    """Iterates over the station record *stream*, returning an
    iterator giving an annual anomaly series for each station as
    (dict, anoms).  The algorithm is as follows: compute monthly
    averages, then monthly anomalies, then seasonal anomalies (means
    of monthly anomalies for at least two months) then annual
    anomalies (means of seasonal anomalies for at least three
    seasons).
    """

    for (dict, series) in stream:
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
            # Seasons are Dec-Feb, Mar-May, Jun-Aug, Sep-Nov.  Dec from previous year.
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
                v = (10.0 * sum(season_anomalies)) / len(season_anomalies)
                annual_anoms.append(int(round(v)))
                if first is None:
                    first = y
                last = y
            else:
                annual_anoms.append(BAD)
        
        if first is not None:
            dict['first'] = first + dict['begin']
            dict['last'] = last + dict['begin']
            yield (dict, annual_anoms[first: last+1])

def is_rural(dict):
    """Test whether the station described by *dict* is rural.
    """
    return dict['US-brightness'] == '1' or (dict['US-brightness'] == ' ' and dict['pop'] == 'R')

class Struct(object):
    pass

# "Global" parameters of the urban_adjustments phase.
g = Struct()
g.xcrit = 2.0 / 3.0
g.x0 = 1950
g.lshort = 7
g.slplim = 1.0
g.slpx = 0.5

g.log = None

def urban_adjustments(anomaly_stream):
    """Takes an iterator of station annual anomaly records (*dict*,
    *series*) and produces an iterator of urban adjustment parameters
    records.  The urban adjustment parameters describe a two-part
    linear fit to the difference in annual anomalies between an urban
    station and the set of nearby rural stations.

    The algorithm is essentially as follows:

    For each urban station:
        1. Find all the rural stations within a fixed radius;
        2. Combine the annual anomaly series for those rural stations, in
           order of valid-data count;
        3. Calculate a two-part linear fit for the difference between
           the urban annual anomalies and this combined rural annual anomaly;
        4. Yield the parameters of this linear fit.

        If there are not enough rural stations, or the combined rural
        record does not have enough overlap with the urban record, try
        a second time for this urban station, with a larger radius.
        If there is still not enough data, do not produce an
        adjustment record.
     """

    g.log = open('log/PApars.GHCN.CL.1000.20.log','w')

    full_radius = 1000.0
    minimum_overlap = 20
    minimum_rural_neighbors = 3

    g.ncrit = 20

    g.nstat = g.sumch = g.nyrs = g.nstap = g.nyrsp = g.sumchp = 0
    g.nsl1 = g.nsl2 = g.ndsl = g.nswitch = g.nsw0 = g.nok = g.nokx = g.nsta = 0
    g.nshort = 0

    f = [0.0] * 900
    x = [0.0] * 900

    # Open output files for some extra logging
    #
    # Station usage stats
    f77 = open("log/PApars.statn.use.GHCN.CL.1000.20", "w")
    # Isolated urban stations
    f79 = open("log/PApars.noadj.stations.list", "w")

    last_year = int(open('work/GHCN.last_year', 'r').read().strip())
    iyoff = 1879
    iyrm = last_year - iyoff

    rural_stations = []
    urban_stations = []
 
    pi180 = math.pi / 180.0

    for (dict, anomalies) in anomaly_stream:
        length = len(anomalies)
        # convert anomalies back from int to float;
        # TODO: preserve as float between annual_anomalies and here.
        for i in range(length):
            if valid(anomalies[i]):
                anomalies[i] *= 0.1

        # throw away some precision in lat/lon for compatibility with Fortran
        lat = 0.1 * math.floor(dict['lat']* 10 + 0.5)
        lon = 0.1 * math.floor(dict['lon']* 10 + 0.5)

        d = Struct()
        d.anomalies = anomalies
        d.cslat = math.cos(lat * pi180)
        d.snlat = math.sin(lat * pi180)
        d.cslon = math.cos(lon * pi180)
        d.snlon = math.sin(lon * pi180)
        d.id = dict['id']
        d.first_year = dict['first'] - iyoff
        d.last_year = d.first_year + length - 1
        d.uses = 0
        if is_rural(dict):
            rural_stations.append(d)
        else:
            urban_stations.append(d)

    nstau = len(urban_stations)        # total number of bright / urban or dim / sm.town stations
    nstar = len(rural_stations)        # total number of dark / rural stations
    g.log.write(" number of rural/urban stations %11d %11d\n" % (nstar, nstau))

    # Sort the rural stations according to the length of the time record
    # (ignoring gaps).
    for i, station in enumerate(rural_stations):
        station.recLen = len([v for v in station.anomalies if valid(v)])
        station.index = i
        g.log.write(" rural station: %11d  id: %s  #ok %11d\n" % (
                    i + 1, station.id, station.recLen))
    rural_stations.sort(key=lambda s:s.recLen)
    rural_stations.reverse()
    for i, station in enumerate(rural_stations):
        g.log.write(" rural station: %11d  id: %s  #ok %11d\n" % (
                    i + 1, station.id, station.recLen))

    # Combine time series for rural stations around each urban station
    for us in urban_stations:
        iyu1 = us.first_year + iyoff - 1           # subtract 1 for a possible partial yr
        iyu2 = us.last_year + iyoff + 1            # add 1 for partial year

        usingFullRadius = False
        dropStation = False
        needNewNeighbours = True
        while True:
            if needNewNeighbours:
                if usingFullRadius:
                    radius = full_radius
                else:
                    radius = full_radius / 2
                neighbors = get_neighbours(us, rural_stations, radius)
                if not neighbors:
                    if usingFullRadius:
                        dropStation = True
                        g.log.write(' no rural neighbors for %s\n' % us.id)
                        f79.write(" no rural neighbors for %s\n" % (us.id))
                        break
                    usingFullRadius = True
                    needNewNeighbours = True
                    continue

                counts, urban_series, combined = combine_neighbors(us, iyrm, iyoff, neighbors, minimum_overlap)
                iy1 = 1
                needNewNeighbours = False

            quorate_count, first, last, quorate_period_total, length = prepare_series(
                iy1, iyrm, combined, urban_series, counts, f, iyoff, x, minimum_rural_neighbors)

            if quorate_count < g.ncrit:
                if usingFullRadius:
                    f79.write("%s  good years: %4d   total years: %4d"
                              " too little rural-neighbors-overlap"
                              " - drop station 9999\n" % (
                        us.id, quorate_count, last - first + 1))
                    dropStation = True
                    break
                usingFullRadius = True
                needNewNeighbours = True
                continue

            if float(quorate_count) >= g.xcrit * (last - first + 1. - .1):
                break

            # not enough good years for the given range (<66%)
            # the  - 0.1 is to prevent equality with potential uncertainty
            # due to hardware rounding or compiler behaviour.
            # nick barnes, ravenbrook limited, 2008 - 09 - 10
            # try to save cases in which the gaps are in the early part:
            iy1 = int(last - (quorate_count - 1) / g.xcrit)
            if iy1 < first + 1:
                iy1 = first + 1                  # avoid infinite loop
            f79.write("%s drop early years %4d-%4d\n" % (
                    us.id, 1 + iyoff, iy1 - 1 + iyoff))

        if dropStation:
            continue

        #===  c subtract urban station and call a curve fitting program
        fit = getfit(length, x, f)
        # find extended range
        iyxtnd = int(round(quorate_count / g.xcrit)) - (last - first + 1)
        g.log.write(" possible range increase %11d %11d %11d\n" % (
                    iyxtnd, quorate_count, last - first + 1))
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

        flag = flags(fit, first + iyoff, last + iyoff)
        yield (us.id, fit, first + iyoff, last + iyoff, n1x, n2x, flag)

    nuse = 0
    for rs in rural_stations:
        if rs.uses > 0:
            f77.write(" used station  %s %11d  times\n" % (
                rs.id, rs.uses))
            nuse += 1
    f77.write("%12d  rural stations were used\n" % (nuse))
    g.log.write(" %-10s %4d %10.7f     %10.7f\n" % (
                "all", g.nsta,-g.sumch/g.nsta,-10.*g.sumch/g.nyrs))
    g.log.write(" %-11s %8d  %10.7f    %10.7f\n" % (
                "urb warm", g.nstap,-g.sumchp/g.nstap,-10.*g.sumchp/g.nyrsp))
    g.log.write(" %-11s %11d %11d %11d %11d %11d %11d\n" % (
                "# short,sl1,sl2,dsl,ok", g.nshort,g.nsl1,g.nsl2,g.ndsl,g.nok,g.nokx))
    g.log.write(" %-11s  %11d %11d\n" % ("switches: all , else ok", g.nswitch,g.nsw0))



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


def combine_neighbors(us, iyrm, iyoff, neighbors, minimum_overlap):
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
    urban_series = [BAD] * iyrm
    combined = [BAD] * iyrm

    urban_series[us.first_year - 1:us.last_year] = us.anomalies
    g.log.write("urb stnID:%s # rur:%4d ranges:%5d%5d.\n" % (
                us.id, len(neighbors), us.first_year + iyoff, us.last_year + iyoff))

    # start with the neighbor with the longest time record
    rs = neighbors[0]
    rs.uses += 1

    combined[rs.first_year - 1:rs.last_year] = rs.anomalies
    for m in range(len(rs.anomalies)):
        if valid(rs.anomalies[m]):
            weights[m + rs.first_year - 1] = rs.weight
            counts[m + rs.first_year - 1] = 1
    g.log.write("longest rur range:%5d-%4d%6d%s\n" % (
                rs.first_year + iyoff, rs.last_year + iyoff, rs.recLen, rs.id))

    # add in the remaining stations
    for i, rs in enumerate(neighbors[1:]):
        g.log.write("add stn%5d range:%5d-%4d %5d %s\n" % (
                    i + 2, rs.first_year + iyoff, rs.last_year + iyoff, rs.recLen,
                    rs.id))
        dnew = [BAD] * iyrm
        dnew[rs.first_year - 1: rs.last_year] = rs.anomalies
        nsm, ncom = cmbine(combined, weights, counts, dnew, rs.first_year, rs.last_year, rs.weight, minimum_overlap)
        g.log.write(" data added:  %11d  overlap: %11d  years\n" % (nsm, ncom))
        if nsm != 0:
            rs.uses += 1

    return counts, urban_series, combined


def prepare_series(iy1, iyrm, combined, urban_series, counts, f, iyoff, x, minimum_rural_neighbors):
    """Prepares for the linearity fitting by populating the arrays *f*
    and *x* with coordinates.  *x* gets a year number (based at
    *g.x0*), and *f* gets the difference between the combined rural
    station anomaly series *combined* and the urban station series
    *urban_series*.  The arrays are only populated with valid years,
    from the first quorate year onwards.  A valid year is one in which
    both the urban station and the combined rural series have valid
    data.  A quorate year is a year in which there are at least
    *minimum_rural_neighbors* contributing (obtained from the *counts*
    series).

    Also returns a 5-tuple: (*c*, *f*, *l*, *t*, *v*). *c* is a count
    of the valid quorate years.  *f* is the first such year.  *l* is
    the last such year.  *t* is the total of the *f* series between
    those two years inclusive.  *v* is the number of entries in the
    *f* and *x* series between those two years inclusive.
    """
    total = first = last = i = quorate_count = length = quorate_period_total = 0

    for iy in xrange(iy1 - 1, iyrm):
        if valid(combined[iy]) and valid(urban_series[iy]):
            if counts[iy] >= minimum_rural_neighbors:
                last = iy + 1
                quorate_count += 1
                if first == 0:
                    first = iy + 1
            if quorate_count <= 0:
                continue

            f[i] = combined[iy] - urban_series[iy]
            total += f[i]
            x[i] = iy + iyoff + 1 - g.x0
            i += 1
            if counts[iy] >= minimum_rural_neighbors:
                 length = i
                 quorate_period_total = total

    return quorate_count, first, last, quorate_period_total, length

    
def cmbine(combined, weights, counts, data, first, last, weight, minimum_overlap):
    """Adds the array *data* with weight *weight* into the array of
    weighted averages *combined*, with total weights *weights* and
    combined counts *counts* (that is, entry *combined[i]* is the
    result of combining *counts[i]* values with total weights
    *weights[i]*).  Adds the computed bias between *combined* and
    *data* before combining.

    Only combines in the range [*first*, *last*); only combines valid
    values from *data*, and if there are fewer than *minimum_overlap*
    entries valid in both arrays then it doesn't combine at all.

    Returns (*n*,*c*), where *n* is the number of entries combined and
    *c* was the number of values valid in both arrays.

    Note: if *data[i]* is valid and *combined[i]* is not, the weighted
    average code runs and still produces the right answer, because
    *weights[i]* will be zero.
    """
    nsm = sumn = ncom = 0
    avg_sum = 0.0
    a, b = first - 1, last
    for v_avg, v_new in itertools.izip(combined[a:b], data[a:b]):
        if invalid(v_avg) or invalid(v_new):
            continue
        ncom = ncom + 1
        avg_sum += v_avg
        sumn += v_new

    if ncom < minimum_overlap:
        return nsm, ncom
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
        nsm += 1

    return nsm, ncom


def getfit(length, x, f):
    """ Finds the best two-part linear fit between *x* and *f*, and
    returns the fit parameters.
    """

    # Todo: incorporate trend2 into this.
    nhalf = length / 2
    rmsmin = 1.e20

    for n in xrange(6, length - 5 + 1):
        xknee = x[n - 1]
        sl1, sl2, rms, sl = trend2(x, f, length, xknee, 2)
        
        if rms < rmsmin:
             rmsmin = rms
             xmin = xknee + g.x0
             fit = (sl1, sl2, xmin, sl)

    return fit


def trend2(xc, a, dataLen, xmid, min):
    """Finds a fit to the data *xc[]*, *a[]*, using regression
    analysis, by a line with a change in slope at *xmid*. Returned is
    a 4-tuple (*sl1*, *sl2*, *rms*, *sl*): the left-hand slope, the
    right-hand slope, the RMS error, and the slope of an overall linear fit.
    """

    # Todo: incorporate into getfit.
    count0 = count1 = 0
    sx0 = sx1 = 0
    sxx0 = sxx1 = 0
    sxa0 = sxa1 = 0

    sa = 0.0
    saa = 0.0

    for n in xrange(dataLen):
        if invalid(a[n]):
            continue
        x = xc[n] - xmid
        v_a = a[n]
        sa += v_a
        saa += v_a ** 2
        if x > 0.0:
            count1 += 1
            sx1 += x
            sxx1 += x ** 2
            sxa1 += x * v_a
        else:
            count0 += 1
            sx0 += x
            sxx0 += x ** 2
            sxa0 += x * v_a

    if count0 < min or count1 < min: 
       return BAD, BAD, BAD, BAD

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

def flags(fit, iy1, iy2):
    """Calculates flags concerning a two-part linear fit.  Also
    accumulate statistics.  Not well documented or understood, but
    note that in adj(), below, the two-part fit will be disregarded if 
    the flag value is not either 0 or 100.
    """

    (sl1, sl2, knee, sl) = fit
    g.nsta += 1
    g.sumch += sl * (iy2 - iy1 + 1)
    g.nyrs += (iy2 - iy1 + 1)
    if sl < 0.0:
        g.nstap += 1
        g.nyrsp += (iy2 - iy1 + 1)
        g.sumchp += sl * (iy2 - iy1 + 1)

    # classify : iflag: +1 for short legs etc
    iflag = 0
    if knee < iy1 + g.lshort or knee > iy2 - g.lshort:
        iflag += 1
        g.nshort += 1
    if abs(sl1) > g.slplim:
        iflag += 20
        g.nsl1 += 1
    if abs(sl2) > g.slplim:
        iflag += 10
        g.nsl2 += 1
    if abs(sl2 - sl1) > g.slplim:
        iflag += 100
        g.ndsl += 1
    if abs(sl2 - sl1) > g.slpx:
        iflag += 100

    if iflag == 0:
        g.nok += 1
    if iflag == 100:
        g.nokx += 1
    if sl1 * sl2 < 0.0 and abs(sl1) > 0.2 and abs(sl2) > 0.2:
        iflag += 1000
        g.nswitch += 1
    if iflag == 1000:
        g.nsw0 += 1
    return iflag

def apply_adjustments(stream, adjustments):
    """Applies the urban adjustment records from the iterator
    *adjustments* to the station records in *stream*.  Returns an
    iterator of adjusted station records.

    Rural stations are passed unchanged.  Urban stations without an
    adjustment record are discarded.  Urban stations with an
    adjustment record are adjusted accordingly, to remove the modelled
    urban-heat-island effect.

    An urban adjustment record describes linear and two-part linear
    fits to the difference between an urban station annual anomaly
    series and the combination of annual anomaly series for nearby
    rural stations.  The linear fit is to allow for a linear urban
    heat-island effect at the urban station.  The two-part linear fit
    is to allow for a model of urbanization starting at some point
    during the time series.
    """

    log = open('log/padjust.log', 'w')
    first_year = 1880
    adj_dict = {}
    for adjustment in adjustments:
        id = adjustment[0]
        adj_dict[id] = adjustment

    for (dict, series) in stream:
        report_name = "%s%c%c%c%3s" % (
                      dict['name'], dict['US-brightness'],
                      dict['pop'], dict['GHCN-brightness'], dict['id'][:3])
        report_station = "station   %9d" % int(dict['id'][3:])
        if adj_dict.has_key(dict['id']):
            # adjust
            m1 = dict['min_month'] - first_year * 12 + 1  # number of first valid month
            m2 = dict['max_month'] - first_year * 12 + 1  # number of last valid month
            offset = dict['min_month'] - dict['begin'] * 12 # index of first valid month
            (_, fit, iy1, iy2, iy1e, iy2e, flag) = adj_dict[dict['id']]
            a, b = adj(first_year, dict, series, fit, iy1e, iy2e, iy1, iy2, flag, m1, m2, offset)
            # a and b are numbers of new first and last valid months
            log.write(" %s  %s saved\n" % (report_station, report_name))
            log.write(" %s  %s adjusted %s %s\n" % (report_station, report_name, (m1, m2), (a, b)))
            aa = a - m1
            bb = b - a + 1
            log.write("ADJ %s %s %s %s %s\n" % (report_name, len(series), len(series), aa, bb))
            series = series[aa + offset:aa + offset + bb]
            dict['begin'] = ((a-1) / 12) + first_year
            dict['first'] = dict['begin']
            dict['min_month'] = a-1 + first_year * 12
            dict['max_month'] = b-1 + first_year * 12
            dict['end'] = ((b-1) / 12) + first_year
            dict['last'] = dict['end']
            log.write("%s %s\n" % (len(series), len(series)))
            yield(dict, series)
            log.write(" %s  %s saved\n" % (report_station, report_name))
        else:
            if is_rural(dict):
                # rural station: not adjusted
                offset = dict['min_month'] - dict['begin'] * 12 # index of first valid month
                length = dict['max_month'] - dict['min_month'] + 1
                series = series[offset : offset + length]
                dict['begin'] = dict['first']
                dict['end'] = dict['last']
                yield (dict, series)
                log.write(" %s  %s saved\n" % (report_station, report_name))
            else:
                # unadjusted urban station: skip
                log.write(" %s  %s skipped\n" % (report_station, report_name))

def adj(first_year, dict, series, fit, iy1, iy2, iy1a, iy2a, iflag, m1, m2, offset):
    (sl1, sl2, knee, sl0) = fit
    if iflag not in (0, 100):
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
        iadj = int(round((iya - knee) * sl - (iy2a - knee) * sl2))
        for m in range(m0, m0 + 12):
            mIdx = m - base
            if mIdx < 0:
                continue
            if m >= m1o and m <= m2o and valid(series[mIdx + offset]):
                if m1 < 0:
                    m1 = m
                series[mIdx+offset] = series[mIdx+offset] + iadj
                m2 = m

        m0 = m0 + 12

    return m1, m2

def results_to_text(input_name, output_name):
    """ Useful for debugging.  Takes a trimmed data file
    *input_name* and produces a text equivalent called *output_name*.
    """

    infile = fort.open(input_name, "rb")
    outfile = open(output_name, 'w')

    header = ccc_binary.CCHeader(infile.readline())
    outfile.write("'%s' %s\n" % (header.title, header.info))

    m1 = header.info[0]
    m2 = header.info[8]

    for s in infile:
        nMonths = m2 - m1 + 1
        rec = ccc_binary.CCRecord(nMonths, data=s)
        outfile.write("id %d lat %d lon %d height %d name %s m1 %d m2 %d\n" % (
                rec.ID, rec.Lat, rec.Lon, rec.iht, rec.name, rec.m1, rec.m2))
        outfile.write("%s\n" % rec.idata)
        m1 = rec.m1
        m2 = rec.m2

def write_results(stream, filename):
    """ Takes a *stream* of step2 result records and writes them to
    a Fortran-format binary file named *filename*.
    """

    first_year = 1880
    last_year = int(open('work/GHCN.last_year', 'r').read().strip())
    MTOT = 12 * (last_year - first_year + 1)
    header = ccc_binary.CCHeader()
    header.title = 'GHCN V2 Temperatures (.1 C)'
    header.info[1] = 1
    header.info[2] = 6
    header.info[3] = MTOT
    header.info[4] = MTOT+15
    header.info[5] = 1880
    header.info[6] = 9999
    header.info[7] = -9999

    outfile = fort.open(filename, "wb")

    header_written = False

    for (dict, series) in stream:
        m1 = dict['min_month'] - first_year * 12 + 1
        m2 = dict['max_month'] - first_year * 12 + 1
        if header_written:
            rec.m1 = m1
            rec.m2 = m2
            outfile.writeline(rec.binary)
        else:
            header.info[0] = m1
            header.info[8] = m2
            outfile.writeline(header.binary)
            header_written = True
        rec = ccc_binary.CCRecord(m2 - m1 + 1)
        rec.idata = series
        rec.Lat = math.floor(dict['lat'] * 10 + 0.5)
        rec.Lon = math.floor(dict['lon'] * 10 + 0.5)
        rec.iht = int(dict['elevs'])
        rec.ID = int(dict['id'][3:])
        rec.name = "%s%c%c%c%3s" % (dict['name'], dict['US-brightness'],
                                    dict['pop'], dict['GHCN-brightness'], dict['id'][:3])

    rec.m1 = 9999
    rec.m2 = 9999
    outfile.writeline(rec.binary)
    
def main(argv=None):
    if argv is None:
        argv = sys.argv

    data = invnt('Ts.GHCN.CL', read_text())
    # At this point we may need to reorder the data so that all the
    # stations between +60.1 and +90.0 comes first, then all the
    # stations between +30.1 and +60.0 come next, and so on.  Thus
    # reflecting how GISTEMP re-orders them when they are split into 6
    # files.  Doing the reordering makes a tiny amount of difference,
    # see http://code.google.com/p/ccc-gistemp/issues/detail?id=25
    #
    # But if you feel like doing this, you'll need to look at the, now
    # deleted, split_binary.py program to see exactly how the split
    # happens.

    (data_1, data_2) = itertools.tee(data, 2)

    anomalies = annual_anomalies(data_1)
    adjustments = urban_adjustments(anomalies)
    
    adjusted = invnt('Ts.GHCN.CL.PA', apply_adjustments(data_2, adjustments))
    write_results(adjusted, 'work/Ts.GHCN.CL.PA')

if __name__ == '__main__':
    main()

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
    """Test for invalid datum (equal to the BAD value).
    """

    return x == BAD

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

def toANNanom(stream):
    """ Iterates over the station record *stream*, returning an
    iterator giving an annual anomaly series for each station as
    (dict, anoms).  The dict is passed through unchanged from the
    station record.  The algorithm is as follows: compute monthly
    averages, then monthly anomalies, then seasonal anomalies - means
    of monthly anomalies for at least two months - then annual
    anomalies - means of seasonal anomalies for at least three
    seasons.
    """

    for (dict, series) in stream:
        av = []
        for m in range(12):
            month_data = filter(valid, series[m::12])
            # neglect December of final year, as we do not use its season.
            if m == 11 and valid(series[-1]):
                month_data = month_data[:-1]
            av.append(float(sum(month_data)) / len(month_data))
        annual_anomalies = []
        # Create groups of seasonal deviations from the month averages.
        first = None
        for y in range(len(series)/12):
            total = [0.0] * 4
            count = [0] * 4
            # Could probably be smarter here by using range(-1,11)
            for m in range(12):
                if m == 11: # Take December value from the previous year
                    year_index = y-1
                else:
                    year_index = y
                if year_index >= 0:
                    datum = series[year_index*12 + m]
                    if valid(datum):
                        season = (m+1) // 3 # season number 0-3
                        if season == 4:
                            season = 0
                        total[season] += datum - av[m]
                        count[season] += 1
            season_anomalies = []
            for s in range(4):
                if count[s] > 1:
                    season_anomalies.append(total[s]/count[s])
            if len(season_anomalies) > 2:
                v = (10.0 * sum(season_anomalies)) / len(season_anomalies)
                annual_anomalies.append(int(round(v)))
                if first is None:
                    first = y
                last = y
            else:
                annual_anomalies.append(BAD)
        
        if first is not None:
            dict['first'] = first + dict['begin']
            dict['last'] = last + dict['begin']
            yield (dict, annual_anomalies[first: last+1])

# PApars.

def is_rural(dict):
    return (dict['US-brightness'] == ' ' and dict['pop'] == 'R') or dict['US-brightness'] == '1'

class Struct(object):
    pass

# Global data
g = Struct()

# Input parameters (# of input files, time period)
g.iyrm0 = 3000 - 1880 + 1

# lower overlap limit, min.coverage of approx.range
g.ncrit = 20
g.nrurm = 3
g.xcrit = 2.0 / 3.0

# Work array sizes
g.nstam = 8000
g.msize = 700000

g.pi180 = math.pi / 180.0
g.x0 = 1950

g.rngbrf = 0.0
g.nlap = 0
g.rngbrh = 0.0
g.rbyrcf = 0.0
g.rbyrch = 0.0
g.cscrif = 0.0
g.cscrih = 0.0
g.mbad = 0
g.xbad = 0.0
g.rdata = []
g.nuseid = []

g.lshort = 7
g.slplim = 1.0
g.slpx = 0.5

g.log = None

def PApars(rngbrf, nlap, anomaly_stream):
    g.log = open('log/PApars.GHCN.CL.1000.20.log','w')
    g.rngbrf = rngbrf
    g.nlap = nlap

    g.nstat = g.sumch = g.nyrs = g.nstap = g.nyrsp = g.sumchp = 0
    g.nsl1 = g.nsl2 = g.ndsl = g.nswitch = g.nsw0 = g.nok = g.nokx = g.nsta = 0
    g.nshort = 0

    f = [0.0] * 900
    x = [0.0] * 900
    w = [0.0] * 900

    g.rngbrh = g.rngbrf / 2
    g.rbyrcf = earth.radius / g.rngbrf
    g.rbyrch = earth.radius / g.rngbrh
    g.cscrif = math.cos(g.rngbrf / earth.radius)
    g.cscrih = math.cos(g.rngbrh / earth.radius)
    g.rural_stations = []
    g.urban_stations = []

    # Open output files for some extra logging
    #
    # Combination info
    f66 = open("log/PApars.statn.log.GHCN.CL.1000.20", "w")
    # Station usage stats
    f77 = open("log/PApars.statn.use.GHCN.CL.1000.20", "w")
    # Isolated urban stations
    f79 = open("log/PApars.noadj.stations.list", "w")

    last_year = int(open('work/GHCN.last_year', 'r').read().strip())
    iyoff = 1879
    iyrm = last_year - iyoff
    nyrsin = iyrm
    g.mbad = 9999
    g.xbad = float(9999)
    current_index = 1

    g.rural_stations = []
    g.urban_stations = []
 
    idata = []
    for (dict, anomalies) in anomaly_stream:
        idata.extend(anomalies)
        length = len(anomalies)

        # throw away some precision in lat/lon for compatibility with Fortran
        lat = 0.1 * math.floor(dict['lat']* 10 + 0.5)
        lon = 0.1 * math.floor(dict['lon']* 10 + 0.5)

        d = Struct()
        d.cslat = math.cos(lat * g.pi180)
        d.snlat = math.sin(lat * g.pi180)
        d.cslon = math.cos(lon * g.pi180)
        d.snlon = math.sin(lon * g.pi180)
        d.id = dict['id']
        d.first_index = current_index
        d.last_index = current_index + length - 1
        d.first_month = dict['first'] - iyoff
        if is_rural(dict):
            g.rural_stations.append(d)
        else:
            g.urban_stations.append(d)

        current_index += length

    nstau = len(g.urban_stations)        # total number of bright / urban or dim / sm.town stations
    nstar = len(g.rural_stations)        # total number of dark / rural stations
    g.log.write(" number of rural/urban stations %11d %11d\n" % (nstar, nstau))

    # Convert data to real numbers (ann avgs were multiplied by 10)
    g.rdata = []
    for v in idata:
        if v == g.mbad:
            g.rdata.append(g.xbad)
        else:
            g.rdata.append(0.1 * v)

    # Sort the rural stations according to the length of the time record
    # (ignoring gaps).
    lengths = []
    for i, station in enumerate(g.rural_stations):
        a, b = station.first_index - 1, station.last_index
        station.recLen = len([v for v in idata[a:b] if v != g.mbad])
        station.index = i
    for i, station in enumerate(g.rural_stations):
        g.log.write(" rural station: %11d  id: %s  #ok %11d\n" % (
                    i + 1, station.id, station.recLen))
    g.rural_stations.sort(key=lambda s:s.recLen)
    g.rural_stations.reverse()
    for i, station in enumerate(g.rural_stations):
        g.log.write(" rural station: %11d  id: %s  #ok %11d\n" % (
                    i + 1, station.id, station.recLen))

    # Combine time series for rural stations around each urban station
    g.nuseid = [0] * nstar
    for us in g.urban_stations:
        usingFullRadius = False
        dropStation = False
        needNewNeighbours = True
        while True:
            if needNewNeighbours:
                data = getNeighbours(us, iyoff, full=usingFullRadius)
                is0, combined, rngbr,    iyu1, iyu2, rbyrc = data
                if is0 == 0:
                    if usingFullRadius:
                        dropStation = True
                        g.log.write(' no rural neighbors for %s\n' % us.id)
                        f79.write(" no rural neighbors for %s\n" % (us.id))
                        break
                    usingFullRadius = True
                    needNewNeighbours = True
                    continue

                wt, iwt, urb, avg = func2(us, iyrm, is0, iyoff, rngbr, combined)
                iy1 = 1
                needNewNeighbours = False

            if iy1 == 1:
                f66.write("year dTs-urban dTs-rural StnID=%s\n" % (
                    us.id))
            tmean, n3l, nxy, n3, n3f, nxy3, tm3 = func3(
                iy1, iyrm, avg, urb, iwt, f, iyoff, x, w, f66)

            if n3 < g.ncrit:
                if usingFullRadius:
                    f79.write("%s  good years: %4d   total years: %4d"
                              " too little rural-neighbors-overlap"
                              " - drop station 9999\n" % (
                        us.id, n3, n3l - n3f + 1))
                    dropStation = True
                    break
                usingFullRadius = True
                needNewNeighbours = True
                continue

            if float(n3) >= g.xcrit * (n3l - n3f + 1. - .1):
                break

            # not enough good years for the given range (<66%)
            # the  - 0.1 is to prevent equality with potential uncertainty
            # due to hardware rounding or compiler behaviour.
            # nick barnes, ravenbrook limited, 2008 - 09 - 10
            # try to save cases in which the gaps are in the early part:
            iy1 = int(n3l - (n3 - 1) / g.xcrit)
            if iy1 < n3f + 1:
                iy1 = n3f + 1                  # avoid infinite loop
            f79.write("%s drop early years %4d-%4d\n" % (
                    us.id, 1 + iyoff, iy1 - 1 + iyoff))

        if dropStation:
            continue

        #===  c subtract urban station and call a curve fitting program
        tmean = tm3 / nxy3
        nxy = nxy3
        fit = getfit(nxy, x, f)
        # find extended range
        iyxtnd = int(round(n3 / g.xcrit)) - (n3l - n3f + 1)
        g.log.write(" possible range increase %11d %11d %11d\n" % (
                    iyxtnd, n3, n3l - n3f + 1))
        n1x = n3f + iyoff
        n2x = n3l + iyoff
        if iyxtnd < 0:
            sys.exit('impossible')
        if iyxtnd > 0:
            lxend = iyu2 - (n3l + iyoff)
            if iyxtnd <= lxend:
                 n2x = n2x + lxend
            else:
                 n1x = n1x - (iyxtnd - lxend)
                 if n1x < iyu1:
                     n1x = iyu1
                 n2x = iyu2

        flag = flags(fit, n3f + iyoff, n3l + iyoff)
        yield (us.id, fit, n3f + iyoff, n3l + iyoff, n1x, n2x, flag)

    nuse = 0
    tempRs = sorted(g.rural_stations, key=lambda s: s.index)
    for (used, tempRs) in itertools.izip(g.nuseid, tempRs):
        if used > 0:
            f77.write(" used station  %s %11d  times\n" % (
                tempRs.id, used))
            nuse += 1
    f77.write("%12d  rural stations were used\n" % (nuse))
    log_stats()


def getNeighbours(us, iyoff, full=False):
    combined = []
    if full:
        cscrit, rbyrc, rngbr = g.cscrif, g.rbyrcf, g.rngbrf
        g.log.write(" trying full radius %s\n" % fort.formatFloat(rngbr))
    else:
        cscrit, rbyrc, rngbr = g.cscrih, g.rbyrch, g.rngbrh

    is0 = 0
    for rs in g.rural_stations:
        iyu1 = us.first_month + iyoff - 1           # subtract 1 for a possible partial yr
        iyu2 = iyu1 + us.last_index - us.first_index + 2  # add 1 for partial year

        csdbyr = (rs.snlat * us.snlat + rs.cslat * us.cslat *
                     (rs.cslon * us.cslon  + rs.snlon * us.snlon))

        if csdbyr <= cscrit:
            continue
        dbyrc = 0
        if csdbyr < 1.0:
            dbyrc = rbyrc * math.sqrt(2.0 * (1.0 - csdbyr))
        is0 += 1
        comb = Struct()
        combined.append(comb)
        comb.wti = 1.0 - dbyrc
        comb.isofi = rs.index + 1
        comb.lenis = rs.recLen
        comb.rs = rs

    return is0, combined, rngbr,     iyu1, iyu2, rbyrc


def func2(us, iyrm, is0, iyoff, rngbr, combined):
    # TODO: Should probably make next 4 arrays a struct (at least as a
    #       stepping stone)
    xbad = g.xbad
    wt = [0.0] * iyrm
    iwt = [0] * iyrm
    urb = [xbad] * iyrm
    avg = [xbad] * iyrm
    ioff = us.first_month - us.first_index
    rdata = g.rdata

    if us.last_index + ioff > g.iyrm0:
        sys.exit("stop 231")
    urb[us.first_index - 1 + ioff:us.last_index + ioff] = rdata[us.first_index - 1:us.last_index]
    g.log.write("urb stnID:%s # rur:%4d ranges:%5d%5d%8.0f.\n" % (
                us.id, is0, us.first_month + iyoff, us.last_index + ioff + iyoff, rngbr))

    #****   start with the station with the longest time record
    comb = combined[0]
    rs = comb.rs
    ioff = rs.first_month - rs.first_index
    g.nuseid[comb.isofi - 1] += 1

    if rs.last_index + 1 + ioff > g.iyrm0:
        sys.exit("stop 244")

    wti = comb.wti
    avg[rs.first_index + ioff - 1:rs.last_index + ioff] = rdata[rs.first_index - 1:rs.last_index]
    for m in xrange(rs.first_index - 1, rs.last_index):
        if rdata[m] < xbad:
            wt[m + ioff] = wti
            iwt[m + ioff] = 1
    g.log.write("longest rur range:%5d-%4d%6d%s\n" % (
                rs.first_month + iyoff, rs.last_index + ioff + iyoff, comb.lenis, rs.id))

    #****   add in the remaining stations
    for i, comb in enumerate(combined[1:is0]):
        is_ = comb.isofi
        rs = comb.rs
        ioff = rs.first_month - rs.first_index
        g.log.write("add stn%5d range:%5d-%4d %5d %s\n" % (
                    i + 2, rs.first_month + iyoff, rs.last_index + ioff + iyoff, comb.lenis,
                    rs.id))
        #****       extend the new data into a full series
        dnew = [xbad] * iyrm
        a, b = rs.first_index - 1, rs.last_index
        dnew[a + ioff: b + ioff] = rdata[a:b]
        nf1 = rs.first_month
        nl1 = rs.last_index + ioff
        #****       shift new data, then combine them with current mean
        nsm, ncom = cmbine(avg, wt, iwt, dnew, nf1, nl1, comb.wti)
        g.log.write(" data added:  %11d  overlap: %11d  years\n" % (nsm, ncom))
        if nsm != 0:
            g.nuseid[is_ - 1] += 1

    return wt, iwt, urb, avg


def func3(iy1, iyrm, avg, urb, iwt, f, iyoff, x, w, f66):
    tmean = nxx = n3l = nxy = n3 = n3f = nxy3 = tm3 = 0

    for iy in xrange(iy1 - 1, iyrm):
        if avg[iy] != g.xbad or urb[iy] != g.xbad:
            nxx = nxx + 1
        if not (avg[iy] == g.xbad or urb[iy] == g.xbad):
            if iwt[iy] >= g.nrurm:
                n3l = iy + 1
                n3 = n3 + 1
                if n3f == 0:
                    n3f = iy + 1
            if n3 <= 0:
                continue

            f[nxy] = avg[iy] - urb[iy]
            tmean = tmean + f[nxy]
            x[nxy] = iy + iyoff + 1 - g.x0
            w[nxy] = 1.
            nxy += 1
            if iwt[iy] >= g.nrurm:
                 nxy3 = nxy
                 tm3 = tmean

        if nxx > 0 and iy1 == 1:
            f66.write("%4d %9.2f %9.2f\n" % (iy + iyoff + 1, urb[iy], avg[iy]))

    return tmean, n3l, nxy, n3, n3f, nxy3, tm3

    
def cmbine(avg, wt, iwt, dnew, nf1, nl1, wt1):
    # bias of new data is removed by subtracting the difference
    # over the common domain. then the new data are averaged in.

    # loop over years
    # find means over common domain to compute bias
    nsm = sumn = ncom = 0
    avg_sum = 0.0
    xbad = g.xbad
    a, b = nf1 - 1, nl1
    for v_avg, v_dnew in itertools.izip(avg[a:b], dnew[a:b]):
        if v_avg >= xbad or v_dnew >= xbad:
            continue
        ncom = ncom + 1
        avg_sum += v_avg
        sumn += v_dnew

    if ncom < g.nlap:
        return nsm, ncom
    bias = (avg_sum - sumn) / float(ncom)

    # update period of valid data, averages and weights
    for n in xrange(nf1 - 1, nl1):
        v_dnew = dnew[n]
        if v_dnew >= xbad:
            continue
        wtnew = wt[n] + wt1
        v_wt, wt[n] = wt[n], wtnew
        avg[n] = (v_wt * avg[n] + wt1 * (v_dnew + bias)) / wtnew
        iwt[n] += 1
        nsm += 1

    return nsm, ncom


def getfit(nxy, x, f):
    # Todo: incorporate trend2 into this.
    nhalf = nxy / 2
    rmsmin = 1.e20

    for n in xrange(6, nxy - 5 + 1):
        xknee = x[n - 1]
        sl1, sl2, rms, sl = trend2(x, f, nxy, xknee, 9999., 2)
        
        if rms < rmsmin:
             rmsmin = rms
             xmin = xknee + g.x0
             fit = (sl1, sl2, xmin, sl)

    return fit


def trend2(xc, a, dataLen, xmid, bad, min):
    # finds a fit using regression analysis by a line
    # with a break in slope at Xmid. Returned are the 2 slopes
    # SL1,SL2 provided we have at least MIN data on each side.
    # Linear regression data are also computed.

    # Todo: incorporate into getfit.
    count0 = count1 = 0
    sx0 = sx1 = 0
    sxx0 = sxx1 = 0
    sxa0 = sxa1 = 0

    sa = 0.0
    saa = 0.0

    for n in xrange(dataLen):
        if a[n] == bad:
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
        return bad, bad, bad, bad

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
    if knee < iy1 + g.lshort or knee > iy2 - g.lshort:
        g.nshort += 1
    if abs(sl1) > g.slplim:
        iflag += 20
    if abs(sl1) > g.slplim:
        g.nsl1 += 1
    if abs(sl2) > g.slplim:
        iflag += 10
    if abs(sl2) > g.slplim:
        g.nsl2 += 1
    if abs(sl2 - sl1) > g.slplim:
        iflag += 100
    if abs(sl2 - sl1) > g.slplim:
        g.ndsl += 1
    if abs(sl2 - sl1) > g.slpx:
        iflag += 100

    if iflag == 0:
        g.nok += 1
    if iflag == 100:
        g.nokx += 1
    if sl1 * sl2 < 0.0 and abs(sl1) > 0.2 and abs(sl2) > 0.2:
        iflag += 1000
    if iflag >= 1000:
        g.nswitch += 1
    if iflag == 1000:
        g.nsw0 += 1
    return iflag

def log_stats():
    g.log.write(" %-10s %4d %10.7f     %10.7f\n" % (
                "all", g.nsta,-g.sumch/g.nsta,-10.*g.sumch/g.nyrs))
    g.log.write(" %-11s %8d  %10.7f    %10.7f\n" % (
                "urb warm", g.nstap,-g.sumchp/g.nstap,-10.*g.sumchp/g.nyrsp))
    g.log.write(" %-11s %11d %11d %11d %11d %11d %11d\n" % (
                "# short,sl1,sl2,dsl,ok", g.nshort,g.nsl1,g.nsl2,g.ndsl,g.nok,g.nokx))
    g.log.write(" %-11s  %11d %11d\n" % ("switches: all , else ok", g.nswitch,g.nsw0))


def padjust(stream, adjustments):
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

    miss = 9999
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

    anomalies = toANNanom(data_1)
    adjustments = PApars(1000.0, 20, anomalies)
    
    adjusted = invnt('Ts.GHCN.CL.PA', padjust(data_2, adjustments))
    write_results(adjusted, 'work/Ts.GHCN.CL.PA')

if __name__ == '__main__':
    main()

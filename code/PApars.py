#!/usr/bin/env python
"""Python replacement for code/STEP2/PApars.f

Input files:

    fort.31,... copies/links of ANN.dTs.GHCN.CL.1,...

"""
__docformat__ = "restructuredtext"


import time
import sys
import os
import math
from itertools import izip

import fort
import ccc_binary
import script_support

options = None

# Input parameters (# of input files, time period)
km0 = 1
iyrm0 = km0 * (3000 - 1880 + 1)

# Earth radius, lower overlap limit, min.coverage of approx.range
rearth = 6375.
ncrit = 20
nrurm = 3
xcrit = 2.0 / 3.0

# Work array sizes
nstam = 8000
msize = 700000

pi180 = math.pi / 180.0
x0 = 1950

class Struct(object):
    pass


# Global data - pretends to be Fortran COMMON (to a limited extent).
g = Struct()
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

def timer():
    start = time.time()
    while True:
        yield "%.2f" % (time.time() - start)


def verbose(level, s):
    if options and level <= options.verbose:
        print >>sys.stderr, s


def main(args):
    """The main for this module.

    """
    T = timer()
    verbose(1, "main T=%s" % T.next())

    # The Python version of this program currently requires the two command
    # line arguments. We do not provide default values.
    g.rngbrf = float(args[0])
    g.nlap = int(args[1])

    ts = [0.0] * 900
    f = [0.0] * 900
    x = [0.0] * 900
    w = [0.0] * 900
    yr = [0.0] * 900
    fpar = [0.0] * 20
    rmsp = [0.0] * 20

    g.rngbrh = g.rngbrf / 2
    g.rbyrcf = rearth / g.rngbrf
    g.rbyrch = rearth / g.rngbrh
    g.cscrif = math.cos(g.rngbrf / rearth)
    g.cscrih = math.cos(g.rngbrh / rearth)
    g.isrData = []
    g.isuData = []

    # Open the 6 input files and single output file.
    infiles = [fort.open("work/ANN.dTs.GHCN.CL.%d" % i) for i in range(1, 7)]
    f78 = open("work/PApars.pre-flags", "w")

    # And for some extra logging
    f66 = open("log/PApars.statn.log.GHCN.CL.1000.20", "w")  # combination info
    f77 = open("log/PApars.statn.use.GHCN.CL.1000.20", "w")  # station usage stats
    f79 = open("log/PApars.noadj.stations.list", "w")  # isolated urban stations

    # Read header of the first file.
    header = ccc_binary.CCHeader(data=infiles[0].readline())
    kq = header.info[1]
    if kq != 1:
        sys.exit("program not ready for quantity %s" % kq)

    # The Fortran version supports different inputs, but we support just the
    # one - corresponding to KMO=1 in PApars.f.
    if header.info[2] in [6, 7]:
        sys.exit("input does not have correct value for info(3)")

    ml = header.info[3]
    iyrm = header.info[3]
    nyrsin = header.info[3]
    if iyrm0 < header.info[3]:
        sys.exit("increase iyrm0 to at least %s" % header.info[3])
    iyoff = header.info[5] - 1
    g.mbad = header.info[6]
    last = header.info[6]
    g.xbad = float(g.mbad)

    isu = 0
    isr = 0
    i1snow = 1

    g.isrData = []
    g.isuData = []
 
    idata = []
    for fidx in range(6):
        progress = script_support.countReport(
                50, fmt="file %s: %%d records processed\n" % fidx)
        infiles[fidx].seek(0)
        header = ccc_binary.CCHeader(data=infiles[fidx].readline())
        mf = header.info[0]
        ml = header.info[8]

        recc = 0
        for recordcount, s in enumerate(infiles[fidx]):
            if i1snow + ml - mf > msize:
                sys.exit('error: msize too small')

            nmonths = ml - mf + 1 
            inrec = ccc_binary.CCRecord(nmonths, data=s)
            idata.extend(inrec.idata)
            recc += 1

            mfcur = mf
            length = ml - mf + 1
            mf, ml = inrec.m1, inrec.m2
            lat = 0.1 * inrec.Lat
            lon = 0.1 * inrec.Lon

            # note isr = rural stations, isu = non-rural stations
            if inrec.name[30:32] == " R" or inrec.name[30] == "1":
                d = Struct()
                g.isrData.append(d)
                
                d.mfsr = mfcur
                d.i1sr = i1snow
                d.ilsr = i1snow + length - 1
                d.cslatr = math.cos(lat * pi180)
                d.snlatr = math.sin(lat * pi180)
                d.cslonr = math.cos(lon * pi180)
                d.snlonr = math.sin(lon * pi180)
                d.idr = inrec.ID
            else:
                d = Struct()
                g.isuData.append(d)

                d.mfsu = mfcur
                d.i1su = i1snow
                d.ilsu = i1snow + length - 1
                d.cslatu = math.cos(lat * pi180)
                d.snlatu = math.sin(lat * pi180)
                d.cslonu = math.cos(lon * pi180)
                d.snlonu = math.sin(lon * pi180)
                d.idu = inrec.ID
                d.cc = inrec.name[33:36]

            i1snow = i1snow + length

            if options and options.verbose >= 2:
                progress.next()

            # Short-circuit
            #if len(g.isrData) >=100:
            #    break
        # Short-circuit
        #if len(g.isrData) >= 100:
        #    break

    verbose(1, "Load T=%s" % T.next())

    nstau = len(g.isuData)        # total number of bright / urban or dim / sm.town stations
    nstar = len(g.isrData)        # total number of dark / rural stations
    ldtot = i1snow - 1          # total length of idata used
    print " number of rural/urban stations %11d %11d" % (nstar, nstau)

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
    for i, station in enumerate(g.isrData):
        a, b = station.i1sr - 1, station.ilsr
        station.recLen = len([v for v in idata[a:b] if v != g.mbad])
        station.index = i
    for i, station in enumerate(g.isrData):
        print " rural station: %11d  id: %11d  #ok %11d" % (
                i + 1, station.idr, station.recLen)
    g.isrData.sort(key=lambda s:s.recLen)
    g.isrData.reverse()
    for i, station in enumerate(g.isrData):
        print " rural station: %11d  id: %11d  #ok %11d" % (
                i + 1, station.idr, station.recLen)
    verbose(1, "Sort T=%s" % T.next())

    # Combine time series for rural stations around each urban station
    g.nuseid = [0] * nstar
    for nurb, us in enumerate(g.isuData):
        if nurb % 100 == 0:
            verbose(1, "Combined %s of %s urban stations at T=%s" % (
                nurb, len(g.isuData), T.next()))
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
                        print ' no rural neighbors for %9d' % us.idu
                        f79.write(" no rural neighbors for %9d\n" % (us.idu))
                        break
                    usingFullRadius = True
                    needNewNeighbours = True
                    continue

                wt, iwt, urb, avg = func2(us, iyrm, is0, iyoff, rngbr, combined)
                iy1 = 1
                needNewNeighbours = False

            if iy1 == 1:
                f66.write("year dTs-urban dTs-rural StnID=%9d\n" % (
                    us.idu))
            tmean, n3l, nxy, n3, n3f, nxy3, tm3 = func3(
                iy1, iyrm, avg, urb, iwt, ts, f, iyoff, yr, x, w, f66)

            if n3 < ncrit:
                if usingFullRadius:
                    # print " Drop station", us.cc, " ", us.idu
                    f79.write("%3s%09d  good years: %4d   total years: %4d"
                              " too little rural-neighbors-overlap"
                              " - drop station 9999\n" % (
                        us.cc, us.idu, n3, n3l - n3f + 1))
                    dropStation = True
                    break
                usingFullRadius = True
                needNewNeighbours = True
                continue

            if float(n3) >= xcrit * (n3l - n3f + 1. - .1):
                break

            # not enough good years for the given range (<66%)
            # the  - 0.1 is to prevent equality with potential uncertainty
            # due to hardware rounding or compiler behaviour.
            # nick barnes, ravenbrook limited, 2008 - 09 - 10
            # try to save cases in which the gaps are in the early part:
            iy1 = int(n3l - (n3 - 1) / xcrit)
            if iy1 < n3f + 1:
                iy1 = n3f + 1                  # avoid infinite loop
            # print "%3s%08d drop early years %4d-%4d" % (
            #         us.cc, us.idu, 1 + iyoff, iy1 - 1 + iyoff)
            f79.write("%3s%09d drop early years %4d-%4d\n" % (
                    us.cc, us.idu, 1 + iyoff, iy1 - 1 + iyoff))

        if dropStation:
            continue

        #===  c subtract urban station and call a curve fitting program
        tmean = tm3 / nxy3
        nxy = nxy3
        getfit(nxy, x, f, fpar, rmsp)
        # find extended range
        iyxtnd = int(round(n3 / xcrit)) - (n3l - n3f + 1)
        print " possible range increase %11d %11d %11d" % (
                iyxtnd, n3, n3l - n3f + 1)
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
           
        # write out a table entry for the table of adjustment parameters
        s = ("%3s%09d %8.3f %8.3f %4d %8.3f %8.3f %8.3f %8.3f %8.3f"
                  " %4d-%4d %4d-%4d") % (
                us.cc, us.idu, fpar[0], fpar[1], int(fpar[2] + x0),
                fpar[3], fpar[4], fpar[5],
                rmsp[0], rmsp[1], n3f + iyoff, n3l + iyoff, n1x, n2x)
        f78.write("%s\n" % s)

    verbose(1, "Comb T=%s" % T.next())
    nuse = 0
    tempRs = sorted(g.isrData, key=lambda s: s.index)
    for n, (used, tempRs) in enumerate(izip(g.nuseid, tempRs)):
        if used > 0:
            f77.write(" used station  %11d %11d  times\n" % (
                tempRs.idr, used))
            nuse += 1
    f77.write("%12d  rural stations were used\n" % (nuse))


def getNeighbours(us, iyoff, full=False):
    combined = []
    if full:
        cscrit, rbyrc, rngbr = g.cscrif, g.rbyrcf, g.rngbrf
        print " trying full radius %s" % fort.formatFloat(rngbr)
    else:
        cscrit, rbyrc, rngbr = g.cscrih, g.rbyrch, g.rngbrh

    is0 = 0
    for rs in g.isrData:
        # print " PAO>>> %11d" % (rs.index + 1)
        iyu1 = us.mfsu + iyoff - 1           # subtract 1 for a possible partial yr
        iyu2 = iyu1 + us.ilsu - us.i1su + 2  # add 1 for partial year

        csdbyr = (rs.snlatr * us.snlatu + rs.cslatr * us.cslatu *
                     (rs.cslonr * us.cslonu + rs.snlonr * us.snlonu))

        if csdbyr <= cscrit:
            continue
        # print "COMP>>> %6.4f %6.4f" % (csdbyr, cscrit)
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
    ioff = us.mfsu - us.i1su
    rdata = g.rdata

    if us.ilsu + ioff > iyrm0:
        sys.exit("stop 231")
    urb[us.i1su - 1 + ioff:us.ilsu + ioff] = rdata[us.i1su - 1:us.ilsu]
    print "urb stnID:%9d # rur:%4d ranges:%5d%5d%8.0f." % (
            us.idu, is0, us.mfsu + iyoff, us.ilsu + ioff + iyoff, rngbr)

    #****   start with the station with the longest time record
    comb = combined[0]
    rs = comb.rs
    ioff = rs.mfsr - rs.i1sr
    g.nuseid[comb.isofi - 1] += 1

    if rs.ilsr + 1 + ioff > iyrm0:
        sys.exit("stop 244")

    wti = comb.wti
    avg[rs.i1sr + ioff - 1:rs.ilsr + ioff] = rdata[rs.i1sr - 1:rs.ilsr]
    for m in xrange(rs.i1sr - 1, rs.ilsr):
        if rdata[m] < xbad:
            wt[m + ioff] = wti
        if rdata[m] < xbad:
            iwt[m + ioff] = 1
    print "longest rur range:%5d-%4d%6d%10d" % (
            rs.mfsr + iyoff, rs.ilsr + ioff + iyoff, comb.lenis, rs.idr)

    #****   add in the remaining stations
    for i, comb in enumerate(combined[1:is0]):
        is_ = comb.isofi
        rs = comb.rs
        ioff = rs.mfsr - rs.i1sr
        print "add stn%5d range:%5d-%4d %5d %9d" % (
            i + 2, rs.mfsr + iyoff, rs.ilsr + ioff + iyoff, comb.lenis,
            rs.idr)
        #****       extend the new data into a full series
        dnew = [xbad] * iyrm
        a, b = rs.i1sr - 1, rs.ilsr
        dnew[a + ioff: b + ioff] = rdata[a:b]
        nf1 = rs.mfsr
        nl1 = rs.ilsr + ioff
        #****       shift new data, then combine them with current mean
        nsm, ncom = cmbine(avg, wt, iwt, dnew, nf1, nl1, comb.wti, rs.idr)
        print " data added:  %11d  overlap: %11d  years" % (nsm, ncom)
        if nsm != 0:
            g.nuseid[is_ - 1] += 1

    return wt, iwt, urb, avg


def func3(iy1, iyrm, avg, urb, iwt, ts, f, iyoff, yr, x, w, f66):
    tmean = nxx = n3l = nxy = n3 = n3f = nxy3 = tm3 = 0

    for iy in xrange(iy1 - 1, iyrm):
        if avg[iy] != g.xbad or urb[iy] != g.xbad:
            nxx = nxx + 1
        if not (avg[iy] == g.xbad or urb[iy] == g.xbad):
            if iwt[iy] >= nrurm:
                n3l = iy + 1
                n3 = n3 + 1
                if n3f == 0:
                    n3f = iy + 1
            if n3 <= 0:
                continue

            ts[nxy] = avg[iy] - urb[iy]
            f[nxy] = ts[nxy]
            tmean = tmean + f[nxy]
            yr[nxy] = iy + iyoff + 1
            x[nxy] = yr[nxy] - x0
            w[nxy] = 1.
            nxy += 1
            if iwt[iy] >= nrurm:
                 nxy3 = nxy
                 tm3 = tmean

        if nxx > 0 and iy1 == 1:
            f66.write("%4d %9.2f %9.2f\n" % (iy + iyoff + 1, urb[iy], avg[iy]))

    return tmean, n3l, nxy, n3, n3f, nxy3, tm3

    
def cmbine(avg, wt, iwt, dnew, nf1, nl1, wt1, id):
    # bias of new data is removed by subtracting the difference
    # over the common domain. then the new data are averaged in.

    # loop over years
    # find means over common domain to compute bias
    nsm = sumn = ncom = 0
    avg_sum = 0.0
    xbad = g.xbad
    a, b = nf1 - 1, nl1
    for v_avg, v_dnew in izip(avg[a:b], dnew[a:b]):
        if v_avg >= xbad or v_dnew >= xbad:
            continue
        ncom = ncom + 1
        avg_sum += v_avg
        sumn += v_dnew

    if ncom < g.nlap:
        return nsm, ncom
    bias = (avg_sum - sumn) / float(ncom)

    # update period of valid data, averages and weights
    # print >>sys.stderr, nf1 - 1, nl1
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


def getfit(nxy, x, f, fpar, rmsp):
    nhalf = nxy / 2
    rmsmin = 1.e20

    for n in xrange(6, nxy - 5 + 1):
        xknee = x[n - 1]
        sl1, sl2, yknee, rms, sl, y0, rms0 = trend2(
                x, f, nxy, xknee, 9999., 2, 2)
        
        if rms < rmsmin:
             rmsmin = rms
             xmin = xknee + x0
             fpar[0] = sl1
             fpar[1] = sl2
             fpar[2] = xknee
             fpar[3] = yknee
             fpar[4] = sl
             fpar[5] = y0
             rmsp[0] = rms / nxy
             rmsp[1] = rms0 / nxy


def trend2(xc, a, dataLen, xmid, bad, min1, min2):
    # finds a fit using regression analysis by a line
    # with a break in slope at Xmid. Returned are the 2 slopes
    # SL1,SL2 provided we have at least MIN1,MIN2 data.
    # Linear regression data are also computed (for emergencies)

    kount0 = kount1 = 0
    sx0 = sx1 = 0
    sxx0 = sxx1 = 0
    sxa0 = sxa1 = 0

    sl1 = bad
    sl2 = bad
    ymid = bad
    sa = 0.0
    saa = 0.0

    for n in xrange(dataLen):
        if a[n] == bad:
            continue
        x = xc[n] - xmid
        v_a = a[n]
        sa = sa + v_a
        saa = saa + v_a ** 2
        if x > 0.0:
            kount1 += 1
            sx1 += x
            sxx1 += x ** 2
            sxa1 += x * v_a
        else:
            kount0 += 1
            sx0 += x
            sxx0 += x ** 2
            sxa0 += x * v_a

    ntot = kount0 + kount1
    denom = ntot * sxx0 * sxx1 - sxx0 * sx1 ** 2 - sxx1 * sx0 ** 2
    xnum1 = sx0 * (sx1 * sxa1 - sxx1 * sa) + sxa0 * (ntot * sxx1 - sx1 ** 2)
    xnum2 = sx1 * (sx0 * sxa0 - sxx0 * sa) + sxa1 * (ntot * sxx0 - sx0 ** 2)

    if kount0 < min1 or kount1 < min2:
        return sl1, sl2, ymid, rms, sl, y0, rms0

    sl1 = xnum1 / denom
    sl2 = xnum2 / denom
    ymid = (sa - sl1 * sx0 - sl2 * sx1) / ntot
    rms = ntot * ymid ** 2 + saa - 2 * ymid * (sa - sl1 * sx0 - sl2 * sx1) + sl1 * sl1 * sxx0 + sl2 * sl2 * sxx1 - 2 * sl1 * sxa0 - 2 * sl2 * sxa1

    # linear regression
    sx0 = sx0 + sx1
    sxx0 = sxx0 + sxx1
    sxa0 = sxa0 + sxa1
    sl = (ntot * sxa0 - sa * sx0) / (ntot * sxx0 - sx0 ** 2)
    y0 = (sa - sl * sx0) / ntot
    rms0 = ntot * y0 ** 2 + saa + sl * sl * sxx0 - 2 * y0 * (sa - sl * sx0) - 2 * sl * sxa0

    return sl1, sl2, ymid, rms, sl, y0, rms0


if __name__ == "__main__":
    import optparse
    usage = "usage: %prog [options] RngbrF NLAP"
    parser = script_support.makeParser(usage)
    options, args = script_support.parseArgs(parser, __doc__, (2, 2))
    main(args)


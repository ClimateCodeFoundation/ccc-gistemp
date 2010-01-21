#!/usr/bin/env python
# $URL: https://ccc-gistemp.googlecode.com/svn/trunk/code/invnt.py $
# $Rev: 162 $
#
# step2.py
#
# Clear Climate Code, 2010-01-15

"""Python replacement for code/STEP2/*
"""
__docformat__ = "restructuredtext"


import sys
import os
import array
import math
import itertools

# Clear Climate Code
import earth
import fort
import ccc_binary

iyrbeg = 1880
iyrend = 3001
nrecpy = 12
monm = nrecpy * (iyrend - iyrbeg + 1)

def invnt(fname):
    """This scans Ts binary files and extracts station information.
    The output is written to log/*.station.list
    """

    log = open('log/$s.station.list','w')
    f = fort.open('work/%s' % fname, "rb")
    header = ccc_binary.CCHeader(data=f.readline())
    n1 = header.info[0]
    n2 = header.info[8]

    if header.info[5] < iyrbeg:
        sys.exit("Beginning year too late")
    if header.info[3] > monm:
        sys.exit("Too much data: %s > %s" % (header.info[3], monm))
    mbad = header.info[6]

    if header.info[0] == mbad:
        # There is no data
        return

    for s in f:
        length = n2 - n1 + 1
        rec = ccc_binary.CCRecord(length, data=s)
        log.write("%s %9d %-30s lat,lon (.1deg)%5d%6d %s%s%s cc=%3s\n" % (
                  fname, rec.ID, rec.name[:30], rec.Lat, rec.Lon,
                  rec.name[31], rec.name[30], rec.name[32],
                  rec.name[33:36],
                  ))
        n1 = rec.m1
        n2 = rec.m2
    log.close()


def text_to_binary(lastyr):
    """Python replacement for the code/STEP2/text_to_binary.f

    This processes the input files:

        work/Ts.txt
            Output from the previous stage. This contains temperator records
            for a set of monitoring stations.

        input/v2.inv
            An inventory file, used to provide extra infomation about the
            temperature records in ``work/Ts.txt``.

    WARNING: I may have got some details wrong below.

    The output is ``work/Ts.bin``, where each record contains:

        idata
            One value per month, in units of 0.1 celcius. There a value for
            each month starting from Jan 1880.
        lat, long, ht
            The monitoring station's position as latitutude, longitude and height.
        id
            <TODO: What is the ID?>
        name
            A name for the monitoring station.
        mtot
            The number of months in the record.

    The program also puts writes to the following log files:

        log/short.station.list
            Lists the stations that were dropped from the output due to
            insufficient data.
        - log/station.log
            A summary for all stations. <TODO: Expand and verify>.
    """
    iyear1, lim, multm, ibad = 1880, 20, 2, 9999
    monthly, iok = [0]*12, [0]*12
    lat, lon, id, ht, hto, mmax = 0, 0, [0]*(multm + 1), 0, 0, [0]*(multm + 1)
    name, nameo, line, li, sid, sidi = "", "", "", "", "", ""

    MTOT = 12 * (lastyr - iyear1 + 1)
    idata = [array.array("i", [0] * MTOT) for i in range(multm + 1)]

    header = ccc_binary.CCHeader()
    header.title = 'GHCN V2 Temperatures (.1 C)'
    header.info[0] = 1
    header.info[1] = 1
    header.info[2] = 6
    header.info[3] = MTOT
    header.info[4] = MTOT+15
    header.info[5] = iyear1
    header.info[6] = 9999
    header.info[7] = -9999
    header.info[8] = MTOT

    inPath = "work/Ts.txt"
    try:
        f1 = open(inPath)
    except IOError, exc:
        sys.exit("Could not open %s for reading\n%s" % (inPath, exc))
    f2 = fort.open('work/Ts.bin', "wb")
    f3 = open('input/v2.inv')
    f88 = open('log/short.station.list', "w")
    f99 = open('log/station.log', "w")

    f2.writeline(header.binary)

    ict = 0
    mult = 0
    id1o = -99

    line, sidi = f1.readline(), '           '
    # 10    continue
    rec = ccc_binary.CCRecord(MTOT)
    while line:
        ict = ict + 1
        if ict % 1000 == 0:
            print "%d processed so far" % ict
        lat, lon, sid, ht, name = fort.unpackRecord(line, 2,
                "i4,i5,a12,i4,a36")
        while True:
            if sidi == sid[:11]:
                break
            li = f3.readline()
            if not li:
                sys.stderr.write("Found no inventory entry for\n  %r\n" %
                        line)
                sys.exit("Inventory file %s is too short" % f3.name)
            sidi = li[:11]
        name = name[:30] + li[101] + li[67] + li[100] + li[:3] + name[35:]
        idfull = int(sid[3:12])
        id1 = int(sid[3:11])

        if id1 == id1o:
            mult = mult + 1
            if mult > multm:
                sys.exit("The 'multm' parameter needs to be larger then %d"
                        % multm)
        else:
            if id1o >= 0:
                for m in range(0, mult+1):
                    if mmax[m] >= lim:
                        writeRecord(f2, rec, idata[m], lato, lono, id[m],
                                hto, nameo, MTOT)
                    else:
                        f88.write("%12d  %sdropped\n" % (id[m], nameo))

            mult = 0
            id1o = id1

        id[mult] = idfull
        idata[mult] = array.array("i", itertools.repeat(ibad, MTOT))
        iok = [0] * 12
        mmax[mult] = 0
        monmin = MTOT+1
        monmax = 0

        # F: 20    continue
        while True:
            # ieof = 1
            line = f1.readline()
            if not line or line[0] == ' ':
                # ieof = 0
                break

            iyr, monthly = fort.unpackRecord(line, 1, "i4, 12i5")
            ix = (iyr - iyear1) * 12
            for m in range(12):
                idatum = monthly[m]
                if idatum != 9999:
                    idata[mult][ix + m] = idatum
                    iok[m - 1] = iok[m - 1] + 1
                    mmax[mult] = max(mmax[mult], iok[m - 1])
                    monmin = min(monmin, ix + m + 1)
                    monmax = ix + m + 1

            # F: goto 20
        # F: 30    continue

        f99.write("%12d%12d%12d%12d%12d%12d\n" % (
                mult, id[mult], mmax[mult],lim, monmin, monmax))
        lato = lat
        lono = lon
        hto = ht
        nameo = name

    if id1o >= 0:
        for m in range(0, mult+1):
            if mmax[m] >= lim:
                writeRecord(f2, rec, idata[m], lato, lono, id[m],
                        hto, nameo, MTOT)
            else:
                f88.write("%12d  %sdropped\n" % (id[m], nameo))

    for f in (f1, f2, f3, f88, f99):
        f.close()


def writeRecord(f2, rec, idata, lato, lono, id, hto, nameo, MTOT):
    rec.idata = list(idata)
    rec.Lat = lato
    rec.Lon = lono
    rec.ID = id
    rec.iht = hto
    rec.name = nameo
    rec.m1 = 1
    rec.m2 = MTOT
    f2.writeline(rec.binary)

def trim_binary(input_name, output_name):
    """Trims a binary file, by removing leading and trailing blocks of
    month data, where ``abs(md) > 8000``. The reduced size records are
    written to the output file.

    Note:
        This (currently) tries to closely follow the form of the original
        fortran code in ``code/STEP2/trim_binary.f``. Hence it is too long, too
        nested, etc compared to typical Python code.
    """

    f2 = fort.open(input_name, "rb")
    f3 = ccc_binary.BufferedOutputRecordFile(output_name)

    header = ccc_binary.CCHeader(f2.readline())
    f3.writeRecord(header)

    i4 = header.info[3]
    inRec = ccc_binary.CCRecord(i4)
    marker = header.info[6]

    # TODO: What if m1 == m2 == marker. Should that be:
    #
    # - a fatal condition.
    # - cause an empty record o be written.
    # - cause no record to be written.

    for recordCount, s in enumerate(f2):
        inRec.setFromBinary(s)
        md = inRec.idata[0:i4]
        m1 = marker
        m2 = marker
        for m, v in enumerate(md):
            if abs(v) > 8000:
                md[m] = marker
            else:
                m2 = m  + 1
                if m1 == marker:
                    m1 = m + 1
        rec = f3.lastRecord()
        if recordCount == 0:
            rec.info[0] = m1
            rec.info[8] = m2
        else:
            rec.m1 = m1
            rec.m2 = m2

        # Add a new output record to the buffer.
        outRec = ccc_binary.CCRecord(m2 - m1 + 1)
        outRec.idata = md[m1 - 1: m2]
        outRec.Lat = inRec.Lat
        outRec.Lon = inRec.Lon
        outRec.ID = inRec.ID
        outRec.iht = inRec.iht
        outRec.name = inRec.name
        outRec.m1 = marker
        outRec.m2 = marker
        f3.writeRecord(outRec)

iylast = 3000
i4o = (iylast - 1701 + 1)
i4 = i4o * 12

def yearMonthEnumerate(monthData, m1):
    """Enumerate months and years, with Dec=0 and Nov=11, year=Dec->Nov"""
    yearBase, monthIdx = divmod(m1, 12)
    yearIdx = yearBase
    if m1 % 12 == 0:
        yearBase -= 1
    for d in monthData:
        yield (yearIdx - yearBase, monthIdx), d
        monthIdx += 1
        if monthIdx > 11:
            monthIdx = 0
            yearIdx += 1


def annav(mon, nyrs, iy1, ibad, m1, recordIdx, emuBug=False):
    # Work out the average for each month of the year (Jan, Feb, ..., Dec).
    # First group values for each month of the year.
    if emuBug:
        mon = mon[:-1]
    av = [[] for m in range(12)]
    t = []
    for (yIdx, mIdx), v in yearMonthEnumerate(mon, m1):
        if v != ibad:
            av[mIdx].append(v)
            t.append((yIdx + 1, mIdx + 1, v))
    nYears = yIdx + 1

    # Use the set of value for each month to get an average.
    for m, monthData in enumerate(av):
        if not monthData:
            sys.exit("station too short - impossible, %s: %s" % (m, monthData))
        av[m] = float(sum(monthData)) / len(monthData)

    # Create groups of seasonal deviations from the month averages.
    ss = [[ [], [], [], [] ] for y in range(nYears)]
    for (yIdx, mIdx), v in yearMonthEnumerate(mon, m1):
        if v != ibad:
            sIdx = mIdx // 3
            ss[yIdx][sIdx].append(v - av[mIdx])

    # Average each season's deviation and then average those for each year.
    for year in ss:
        for i, s in enumerate(year):
            year[i] = ibad
            if len(s) > 1:
                year[i] = sum(s) / len(s)

    # Average for the seasons for each year.
    iann = [[] for i in range(nYears)]
    iy1n = iy2n = ibad
    for yIdx, year in enumerate(ss):
        data = [v for v in year if v != ibad]
        iann[yIdx] = ibad
        if len(data) > 2:
            v = (10.0 * sum(data)) / len(data)
            iann[yIdx] = int(round(v))

            iy2n = iy1 + yIdx
            if iy1n == ibad:
                iy1n = iy2n

    return iy1n, iy2n, iann
    

def toANNanom(input_name, output_name, count = None):
    f2 = fort.open(input_name, "rb")
    f3 = ccc_binary.BufferedOutputRecordFile(output_name)
    header = ccc_binary.CCHeader(data=f2.readline())
    ibad = header.info[7 - 1]
    if header.info[1 - 1] == ibad:
        sys.exit("no stations")

    outHdr = header.copy()
    outHdr.info[3 - 1] = 5                         # ann.means (6=mon.means)
    outHdr.info[4 - 1] = header.info[4 - 1] / 12   # length of time series
    outHdr.info[5 - 1] = outHdr.info[4 - 1] + header.info[5 - 1] - header.info[4 - 1]  # length of records
    outHdr.title = "ANNUAL MEAN TEMPERATURE ANOMALIES (.01 C)".ljust(80)

    m1 = header.info[0]
    m2 = header.info[8]

    f3.writeRecord(outHdr)
    doneHeader = False

    for recordCount, s in enumerate(f2):
        nMonths = m2 - m1 + 1
        inRec = ccc_binary.CCRecord(nMonths, data=s)
    
        iy1 = 1 + (m1 - 1) / 12
        mon1 = 12 * (iy1 - 1)   #   =  december of year iy1 - 1
        nyrs = (m2 + 11 - mon1) / 12
        if nyrs - 1 + iy1 > i4o:
            nyrs = i4o + 1 - iy1
        if nyrs < 0:
            sys.exit('no station records  -  impossible')

        i1, i2, iann = annav(inRec.idata[0:nMonths], nyrs, iy1, ibad, m1, recordCount,
                emuBug=(m2 % 12 == 0))
        if count is not None and recordCount > count:
            return
        m1, m2 = inRec.m1, inRec.m2
        if i1 == ibad:
            continue
        if not doneHeader:
            outHdr.info[0] = i1
            outHdr.info[8] = i2
            doneHeader = True
        else:
            rec = f3.lastRecord()
            rec.m1 = i1
            rec.m2 = i2

        nYears = i2 - i1 + 1
        outRec = ccc_binary.CCRecord(nYears)
        outRec.idata = iann[i1-iy1:i2-iy1+1]
        outRec.Lat = inRec.Lat
        outRec.Lon = inRec.Lon
        outRec.ID = inRec.ID
        outRec.iht = inRec.iht
        outRec.name = inRec.name
        outRec.m1 = ibad
        outRec.m2 = ibad
        f3.writeRecord(outRec)


def readRec(f):
    l = f.readline()
    if not l:
        return None, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

    (CCdStationID,slope_l,slope_r,knee,Yknee,slope,Ymid,RMS,RMSl,
            rur_urb,ext_range,flag) = l.split()
    cc, IDc = CCdStationID[:3], int(CCdStationID[3:])
    sl1 = float(slope_l)
    sl2 = float(slope_r)
    knee = int(knee)
    sl0 = float(slope)
    iy1, iy2 = [int(v) for v in rur_urb.split("-")]
    iy1e, iy2e = [int(v) for v in ext_range.split("-")]
    flag = int(flag)
    return cc, IDc, sl1, sl2, knee, sl0, iy1, iy2, iy1e, iy2e, flag

def padjust():
    log = open('log/padjust.log', 'w')
    f1 = open("work/PApars.list")
    header = f1.readline()
    # Read the first adjustment record.
    cc, IDc, sl1, sl2, knee, sl0, iy1, iy2, iy1e, iy2e, flag = readRec(f1)

    inF = fort.open("work/Ts.GHCN.CL")
    outF = fort.open("work/Ts.GHCN.CL.PA", "wb")

    h = ccc_binary.CCHeader(inF.readline())
    ibad = h.info[6]
    if h.info[0] == ibad:
        raise Exception("was break, but drj removed the while loop")
    m1o, m2o = h.info[0], h.info[8]

    # Loop until we find a station we do not want to skip. At that point,
    # we start writing the current output file (i.w. write the header) and
    # adjust the first record if necessary.
    while True:
        count = m2o - m1o + 1
        s = data=inF.readline()
        if not s:
            break
        rec_o = ccc_binary.CCRecord(count, data=s)
        if rec_o.ID != IDc:
            if rec_o.name[30:32] == " R" or rec_o.name[30] == "1":
                outF.writeline(h.binary)
                break
            else: # skip the station
                h.info[0] = rec_o.m1
                h.info[8] = rec_o.m2
                m1o, m2o = rec_o.m1, rec_o.m2
                log.write(" station   %9s  %s skipped\n" % (
                                  rec_o.ID, rec_o.name))
                continue

        else:
            # Record needs adjusting. Do so and get the next adjustment
            # entry.
            a, b = adj(h, rec_o.idata, sl1, sl2, knee, sl0, iy1e, iy2e, iy1, iy2,
                    flag, m1o, m2o)
            log.write( " station   %9s  %s adjusted %s %s\n" % (
                              rec_o.ID, rec_o.name, (m1o, m2o), (a, b)))
            aa = a - m1o
            bb = b - a + 1
            log.write("ADJ-1 %s %s %s %s %s\n" % (rec_o.name, rec_o.count, len(rec_o.idata), aa, bb))
            rec_o.idata[:] = rec_o.idata[aa:aa + bb]
            rec_o.count = bb
            log.write("%s %s\n" % (rec_o.count, len(rec_o.idata)))
            m1o, m2o = a, b
            h.info[0] = m1o
            h.info[8] = m2o
            outF.writeline(h.binary)
            ddd = readRec(f1)
            if ddd[0] is not None:
                cc, IDc, sl1, sl2, knee, sl0, iy1, iy2, iy1e, iy2e, flag = ddd
            break

    # Now read further records.
    while True:
        s = data=inF.readline()
        if not s:
            outF.writeline(rec_o.binary)
            log.write(" station   %9s  %s saved\n" % (rec_o.ID, rec_o.name))
            break

        rec = ccc_binary.CCRecord(rec_o.m2 - rec_o.m1 + 1, data=s)
        if rec.ID != IDc:
            if rec.name[30:32] == " R" or rec.name[30] == "1":
                outF.writeline(rec_o.binary)
                log.write(" station   %9s  %s saved\n" % (rec_o.ID, rec_o.name))
            else:
                rec_o.m1, rec_o.m2 = rec.m1, rec.m2
                log.write(" station   %9s  %s skipped\n" % (rec.ID, rec.name))
                continue

        else:
            a, b = adj(h, rec.idata, sl1, sl2, knee, sl0, iy1e, iy2e, iy1, iy2,
                    flag, rec_o.m1, rec_o.m2)
            log.write(" station   %9s  %s adjusted %s %s\n" % (
                              rec.ID, rec.name, (rec_o.m1, rec_o.m2), (a, b)))
            aa = a - rec_o.m1
            bb = b - a + 1
            log.write("ADJ-2 %s %s %s %s %s\n" %(rec.name, rec.count, len(rec.idata), aa, bb))
            rec.idata[:] = rec.idata[aa:aa + bb]
            rec.count = bb
            log.write("%s %s\n" % (rec.count, len(rec.idata)))
            rec_o.m1, rec_o.m2 = a, b
            outF.writeline(rec_o.binary)
            log.write(" station   %9s  %s saved\n" % (rec_o.ID, rec_o.name))
            IDc = -9999
            ddd = readRec(f1)
            if ddd[0] is not None:
                cc, IDc, sl1, sl2, knee, sl0, iy1, iy2, iy1e, iy2e, flag = ddd

        rec_o = rec.copy()


def adj(h, idata, sl1, sl2, knee, sl0, iy1, iy2, iy1a, iy2a, iflag, m1, m2):
    if iflag not in (0, 100):
        # Use linear approximation
        sl1, sl2 = sl0, sl0

    base = m1

    miss = h.info[6]
    m1o, m2o = m1, m2
    m1 = -100
    m0 = 12 * (iy1 - h.info[5])   # Dec of year iy1
    for iy in xrange(iy1, iy2 + 1):
        sl = sl1
        if iy > knee:
            sl = sl2
        iya = iy
        if iy < iy1a:
            iya = iy1a
        if iy > iy2a:
            iya = iy2a
        iadj = int(round((iya - knee) * sl - (iy2a - knee) * sl2))
        for m in xrange(m0, m0 + 12):
            mIdx = m - base
            if mIdx < 0:
                continue
            if m >= m1o and m <= m2o and idata[mIdx] != miss:
                if m1 < 0:
                    m1 = m
                idata[mIdx] = idata[mIdx] + iadj
                m2 = m

        m0 = m0 + 12

    return m1, m2


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

g.log = None

def PApars(rngbrf, nlap):
    g.log = open('log/PApars.GHCN.CL.1000.20.log','w')
    g.rngbrf = rngbrf
    g.nlap = nlap

    ts = [0.0] * 900
    f = [0.0] * 900
    x = [0.0] * 900
    w = [0.0] * 900
    yr = [0.0] * 900
    fpar = [0.0] * 20
    rmsp = [0.0] * 20

    g.rngbrh = g.rngbrf / 2
    g.rbyrcf = earth.radius / g.rngbrf
    g.rbyrch = earth.radius / g.rngbrh
    g.cscrif = math.cos(g.rngbrf / earth.radius)
    g.cscrih = math.cos(g.rngbrh / earth.radius)
    g.isrData = []
    g.isuData = []

    # Open the input and output files.
    infile = fort.open("work/ANN.dTs.GHCN.CL")
    f78 = open("work/PApars.pre-flags", "w")

    #
    # Open output files for some extra logging
    #
    # Combination info
    f66 = open("log/PApars.statn.log.GHCN.CL.1000.20", "w")
    # Station usage stats
    f77 = open("log/PApars.statn.use.GHCN.CL.1000.20", "w")
    # Isolated urban stations
    f79 = open("log/PApars.noadj.stations.list", "w")

    # Read header of the input file.
    header = ccc_binary.CCHeader(data=infile.readline())
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
    if g.iyrm0 < header.info[3]:
        sys.exit("increase g.iyrm0 to at least %s" % header.info[3])
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
    infile.seek(0)
    header = ccc_binary.CCHeader(data=infile.readline())
    mf = header.info[0]
    ml = header.info[8]

    recc = 0
    for s in infile:
        if i1snow + ml - mf > g.msize:
            sys.exit('error: g.msize too small')

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
            d.cslatr = math.cos(lat * g.pi180)
            d.snlatr = math.sin(lat * g.pi180)
            d.cslonr = math.cos(lon * g.pi180)
            d.snlonr = math.sin(lon * g.pi180)
            d.idr = inrec.ID
        else:
            d = Struct()
            g.isuData.append(d)

            d.mfsu = mfcur
            d.i1su = i1snow
            d.ilsu = i1snow + length - 1
            d.cslatu = math.cos(lat * g.pi180)
            d.snlatu = math.sin(lat * g.pi180)
            d.cslonu = math.cos(lon * g.pi180)
            d.snlonu = math.sin(lon * g.pi180)
            d.idu = inrec.ID
            d.cc = inrec.name[33:36]

        i1snow = i1snow + length

    nstau = len(g.isuData)        # total number of bright / urban or dim / sm.town stations
    nstar = len(g.isrData)        # total number of dark / rural stations
    ldtot = i1snow - 1          # total length of idata used
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
    for i, station in enumerate(g.isrData):
        a, b = station.i1sr - 1, station.ilsr
        station.recLen = len([v for v in idata[a:b] if v != g.mbad])
        station.index = i
    for i, station in enumerate(g.isrData):
        g.log.write(" rural station: %11d  id: %11d  #ok %11d\n" % (
                    i + 1, station.idr, station.recLen))
    g.isrData.sort(key=lambda s:s.recLen)
    g.isrData.reverse()
    for i, station in enumerate(g.isrData):
        g.log.write(" rural station: %11d  id: %11d  #ok %11d\n" % (
                    i + 1, station.idr, station.recLen))

    # Combine time series for rural stations around each urban station
    g.nuseid = [0] * nstar
    for us in g.isuData:
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
                        g.log.write(' no rural neighbors for %9d\n' % us.idu)
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

            if n3 < g.ncrit:
                if usingFullRadius:
                    f79.write("%3s%09d  good years: %4d   total years: %4d"
                              " too little rural-neighbors-overlap"
                              " - drop station 9999\n" % (
                        us.cc, us.idu, n3, n3l - n3f + 1))
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
            f79.write("%3s%09d drop early years %4d-%4d\n" % (
                    us.cc, us.idu, 1 + iyoff, iy1 - 1 + iyoff))

        if dropStation:
            continue

        #===  c subtract urban station and call a curve fitting program
        tmean = tm3 / nxy3
        nxy = nxy3
        getfit(nxy, x, f, fpar, rmsp)
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
           
        # write out a table entry for the table of adjustment parameters
        s = ("%3s%09d %8.3f %8.3f %4d %8.3f %8.3f %8.3f %8.3f %8.3f"
                  " %4d-%4d %4d-%4d") % (
                us.cc, us.idu, fpar[0], fpar[1], int(fpar[2] + g.x0),
                fpar[3], fpar[4], fpar[5],
                rmsp[0], rmsp[1], n3f + iyoff, n3l + iyoff, n1x, n2x)
        f78.write("%s\n" % s)

    nuse = 0
    tempRs = sorted(g.isrData, key=lambda s: s.index)
    for (used, tempRs) in itertools.izip(g.nuseid, tempRs):
        if used > 0:
            f77.write(" used station  %11d %11d  times\n" % (
                tempRs.idr, used))
            nuse += 1
    f77.write("%12d  rural stations were used\n" % (nuse))


def getNeighbours(us, iyoff, full=False):
    combined = []
    if full:
        cscrit, rbyrc, rngbr = g.cscrif, g.rbyrcf, g.rngbrf
        g.log.write(" trying full radius %s\n" % fort.formatFloat(rngbr))
    else:
        cscrit, rbyrc, rngbr = g.cscrih, g.rbyrch, g.rngbrh

    is0 = 0
    for rs in g.isrData:
        iyu1 = us.mfsu + iyoff - 1           # subtract 1 for a possible partial yr
        iyu2 = iyu1 + us.ilsu - us.i1su + 2  # add 1 for partial year

        csdbyr = (rs.snlatr * us.snlatu + rs.cslatr * us.cslatu *
                     (rs.cslonr * us.cslonu + rs.snlonr * us.snlonu))

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
    ioff = us.mfsu - us.i1su
    rdata = g.rdata

    if us.ilsu + ioff > g.iyrm0:
        sys.exit("stop 231")
    urb[us.i1su - 1 + ioff:us.ilsu + ioff] = rdata[us.i1su - 1:us.ilsu]
    g.log.write("urb stnID:%9d # rur:%4d ranges:%5d%5d%8.0f.\n" % (
                us.idu, is0, us.mfsu + iyoff, us.ilsu + ioff + iyoff, rngbr))

    #****   start with the station with the longest time record
    comb = combined[0]
    rs = comb.rs
    ioff = rs.mfsr - rs.i1sr
    g.nuseid[comb.isofi - 1] += 1

    if rs.ilsr + 1 + ioff > g.iyrm0:
        sys.exit("stop 244")

    wti = comb.wti
    avg[rs.i1sr + ioff - 1:rs.ilsr + ioff] = rdata[rs.i1sr - 1:rs.ilsr]
    for m in xrange(rs.i1sr - 1, rs.ilsr):
        if rdata[m] < xbad:
            wt[m + ioff] = wti
        if rdata[m] < xbad:
            iwt[m + ioff] = 1
    g.log.write("longest rur range:%5d-%4d%6d%10d\n" % (
                rs.mfsr + iyoff, rs.ilsr + ioff + iyoff, comb.lenis, rs.idr))

    #****   add in the remaining stations
    for i, comb in enumerate(combined[1:is0]):
        is_ = comb.isofi
        rs = comb.rs
        ioff = rs.mfsr - rs.i1sr
        g.log.write("add stn%5d range:%5d-%4d %5d %9d\n" % (
                    i + 2, rs.mfsr + iyoff, rs.ilsr + ioff + iyoff, comb.lenis,
                    rs.idr))
        #****       extend the new data into a full series
        dnew = [xbad] * iyrm
        a, b = rs.i1sr - 1, rs.ilsr
        dnew[a + ioff: b + ioff] = rdata[a:b]
        nf1 = rs.mfsr
        nl1 = rs.ilsr + ioff
        #****       shift new data, then combine them with current mean
        nsm, ncom = cmbine(avg, wt, iwt, dnew, nf1, nl1, comb.wti, rs.idr)
        g.log.write(" data added:  %11d  overlap: %11d  years\n" % (nsm, ncom))
        if nsm != 0:
            g.nuseid[is_ - 1] += 1

    return wt, iwt, urb, avg


def func3(iy1, iyrm, avg, urb, iwt, ts, f, iyoff, yr, x, w, f66):
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

            ts[nxy] = avg[iy] - urb[iy]
            f[nxy] = ts[nxy]
            tmean = tmean + f[nxy]
            yr[nxy] = iy + iyoff + 1
            x[nxy] = yr[nxy] - g.x0
            w[nxy] = 1.
            nxy += 1
            if iwt[iy] >= g.nrurm:
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


def getfit(nxy, x, f, fpar, rmsp):
    nhalf = nxy / 2
    rmsmin = 1.e20

    for n in xrange(6, nxy - 5 + 1):
        xknee = x[n - 1]
        sl1, sl2, yknee, rms, sl, y0, rms0 = trend2(
                x, f, nxy, xknee, 9999., 2, 2)
        
        if rms < rmsmin:
             rmsmin = rms
             xmin = xknee + g.x0
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

# Constants
g.lshort = 7
g.slplim = 1.0
g.slpx = 0.5


def flags_readRec(f):
    """Read a single record from the input file.

    The input file is text, line oriented, one record per line.

    :Return:
        A tuple of 15 values, as follows:

        cc
            **TBD**
        IDc
            **TBD**
        sl1, sl2
            **TBD**
        knee
            **TBD**
        Yknee
            **TBD**
        sl0, Ymid
            **TBD**
        RMS,RMSl
            **TBD**
        iy1, iy2,
            **TBD**
        iy1e, iy2e
            **TBD**
        l
            **TBD**

    :Param f:
        The open file to read from.

    """
    l = f.readline()
    l = l.rstrip()
    if not l:
        return None, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

    (CCdStationID,slope_l,slope_r,knee,Yknee,slope,Ymid,RMS,RMSl,
            rur_urb,ext_range) = l.split()
    cc, IDc = CCdStationID[:3], int(CCdStationID[3:])
    sl1 = float(slope_l)
    sl2 = float(slope_r)
    knee = int(knee)
    sl0 = float(slope)
    iy1, iy2 = [int(v) for v in rur_urb.split("-")]
    iy1e, iy2e = [int(v) for v in ext_range.split("-")]
    return (cc, IDc, sl1, sl2, knee, Yknee, sl0, Ymid,RMS,RMSl, iy1, iy2,
            iy1e, iy2e, l)


def flags():
    """Python replacement for code/STEP2/flags.f

    Input files:

        work/PApars.pre-flags

    Output files:
        work/PApars.list

    The input file is a text file, with a single record per line. This is created
    by PApars(), which adds the flags value at the end of each record, based on
    the record's contents.

    """ 

    f78 = open("work/PApars.pre-flags")
    f2 = open("work/PApars.list", "w")
    f2.write('CCdStationID  slope-l  slope-r knee    Yknee    slope     Ymid      RMS     RMSl 3-rur+urb ext.range flag\n')
    nstat = sumch = nyrs = nstap = nyrsp = sumchp = 0
    nsl1 = nsl2 = ndsl = nswitch = nsw0 = nok = nokx = nsta = 0
    nshort = 0
    while True:
        data = flags_readRec(f78)
        if data[0] is None:
            break
        cc,id,sl1,sl2,knee,yk,sl,ylin,rms,rms0,iy1,iy2,iy1e,iy2e,line = data
        nsta += 1
        sumch = sumch + sl * (iy2 - iy1 + 1)
        nyrs = nyrs + (iy2 - iy1 + 1)
        if sl < 0.0:
            nstap += 1
            nyrsp += (iy2 - iy1 + 1)
            sumchp += sl * (iy2 - iy1 + 1)

        # classify : iflag: +1 for short legs etc
        iflag = 0
        if knee < iy1 + g.lshort or knee > iy2 - g.lshort:
            iflag += 1
        if knee < iy1 + g.lshort or knee > iy2 - g.lshort:
            nshort += 1
        if abs(sl1) > g.slplim:
            iflag += 20
        if abs(sl1) > g.slplim:
            nsl1 += 1
        if abs(sl2) > g.slplim:
            iflag += 10
        if abs(sl2) > g.slplim:
            nsl2 += 1
        if abs(sl2 - sl1) > g.slplim:
            iflag += 100
        if abs(sl2 - sl1) > g.slplim:
            ndsl += 1

        # TODO: The small offset fixes a Fortran/Python rounding issue, but I
        #       am not happy with it. PAO, Mon Feb 23 2009.
        if abs(sl2 - sl1) + 0.00000001 > g.slpx:
            iflag += 100
       
        if iflag == 0:
            nok += 1
        if iflag == 100:
            nokx += 1
        if sl1 * sl2 < 0.0 and abs(sl1) > 0.2 and abs(sl2) > 0.2:
            iflag += 1000
        if iflag >= 1000:
            nswitch += 1
        if iflag == 1000:
            nsw0 += 1
       
        line += " %4d" % iflag
        f2.write("%s\n" % line)

    g.log.write(" %-10s %4d %10.7f     %10.7f\n" % (
                "all", nsta,-sumch/nsta,-10.*sumch/nyrs))
    g.log.write(" %-11s %8d  %10.7f    %10.7f\n" % (
                "urb warm", nstap,-sumchp/nstap,-10.*sumchp/nyrsp))
    g.log.write(" %-11s %11d %11d %11d %11d %11d %11d\n" % (
                "# short,sl1,sl2,dsl,ok", nshort,nsl1,nsl2,ndsl,nok,nokx))
    g.log.write(" %-11s  %11d %11d\n" % ("switches: all , else ok", nswitch,nsw0))

def main(argv=None):
    if argv is None:
        argv = sys.argv

    year = open('work/GHCN.last_year', 'r').read().strip()
    print "... converting text to binary file, last year = %s" % year
    text_to_binary(int(year))

    # At this point we may need to reorder the Ts.bin/Ts.txt file so
    # that all the stations between +60.1 and +90.0 comes first, then
    # all the stations between +30.1 and +60.0 come next, and so on.
    # Thus reflecting how they get re-ordered when they are split into 6
    # files.
    # Not doing the reordering makes a tiny amount of difference, see
    # http://code.google.com/p/ccc-gistemp/issues/detail?id=25
    # But if you feel like doing, you'll need to look at the, now
    # deleted, split_binary.py program to see exactly how the split
    # happens.

    print "... trimming Ts.bin"
    trim_binary('work/Ts.bin', 'work/Ts.GHCN.CL')

    print "... Making station list"
    invnt('Ts.GHCN.CL')

    print "... Creating annual anomalies"
    toANNanom('work/Ts.GHCN.CL', 'work/ANN.dTs.GHCN.CL')
    PApars(1000.0, 20)
    flags()

    print "... Applying peri-urban adjustment"
    padjust()

    print "... Making station list"
    invnt('Ts.GHCN.CL.PA')

if __name__ == '__main__':
    main()

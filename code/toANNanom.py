#!/usr/bin/env python
"""Python replacement for code/STEP2/toANNanom.f

"""
__docformat__ = "restructuredtext"


import sys
import os

import fort
import ccc_binary
import script_support

options = None

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
    for yIdx, year in enumerate(ss):
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
    

def toANNanom(i):
    fname = "work/Ts.GHCN.CL.%d" % i
    f2 = fort.open(fname, "rb")
    f3 = ccc_binary.BufferedOutputRecordFile("work/ANN.dTs.GHCN.CL.%d" % i)
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

    progress = script_support.countReport(
            50, fmt="File %s: %%d records processed\n" % fname)
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

        if options and options.verbose:
            progress.next()

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


def main(args):
    """The main for this module.

    This simply invokes invntSingleFile for n = 1, 2, ... 6.
    """
    for i in range(1,7):
        toANNanom(i)


if __name__ == "__main__":
    import optparse
    usage = "usage: %prog [options]"
    parser = script_support.makeParser(usage)
    options, args = script_support.parseArgs(parser, __doc__, (0, 0))
    main(args)


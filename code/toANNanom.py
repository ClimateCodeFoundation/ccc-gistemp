#!/usr/bin/env python
"""Python replacement for code/STEP2/toANNanom.f

"""
__docformat__ = "restructuredtext"


import sys
import os

import fort
import ccc_binary
import script_support

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


# Map providing adjustment for a few values that round the wrong way wrt the
# fortran code. The key is (recordIdx, yearIdx, calcValue. The value is an
# adjustment, which should only ever be 1, 0 or -1.
#
# This is a stop-gap measure to make it easier to verify correct operation
# as we develop Python equivalents of each Fortran program.
#
# Note: Entries with an adjusted of zero have, of course, no effect. They
#       exist to provide to make it easy to increase the tolerance used to
#       identify values that may need correcting.
_fixMap = {
    # For output = work/ANN.dTs.GHCN.CL.1
    (240,  9,  -17):   1,
           
    # For output = work/ANN.dTs.GHCN.CL.2
    ( 292,   7,   55):    -1,   # 0x00000037,  frac=0.500000000000
    ( 292,  12,  -70):    -1,   # 0xffffffba,  frac=0.500000000000
    ( 292,  24,   60):    -1,   # 0x0000003c,  frac=0.500000000000
    ( 292,  28,  -20):    -1,   # 0xffffffec,  frac=0.500000000000
    ( 292,  32,  -60):    -1,   # 0xffffffc4,  frac=0.500000000000
    ( 292,  35,  -15):    -1,   # 0xfffffff1,  frac=0.500000000000
    ( 333,   0,    8):     1,   # 0x00000008,  frac=0.500000000000
    ( 558, 104,   85):     0,   # 0x00000055,  frac=0.500012000765
    ( 669,   2,   63):    -1,   # 0x0000003f,  frac=0.500000000000
    ( 716,  11,  -38):     0,   # 0xffffffda,  frac=0.500000000000
    ( 716,  14,   17):     1,   # 0x00000011,  frac=0.500000000000
    ( 716,  22,   97):     1,   # 0x00000061,  frac=0.500000000000
    ( 716,  23,   32):     1,   # 0x00000020,  frac=0.500000000000
    ( 716,  25,  -93):     0,   # 0xffffffa3,  frac=0.500000000000
    ( 716,  34, -198):     0,   # 0xffffff3a,  frac=0.500000000000
    ( 716,  38,   67):     1,   # 0x00000043,  frac=0.500000000000
    ( 716,  39,  132):     1,   # 0x00000084,  frac=0.500000000000
    ( 716,  40,  218):     0,   # 0x000000da,  frac=0.500000000000
    ( 716,  44,  -53):     1,   # 0xffffffcb,  frac=0.500000000000
    ( 749,  12,  -75):    -1,   # 0xffffffb5,  frac=0.500000000000
    ( 749,  16,  -95):    -1,   # 0xffffffa1,  frac=0.500000000000
    (1295,  18,   64):     0,   # 0x00000040,  frac=0.500000000000
    (1295,  22,  140):    -1,   # 0x0000008c,  frac=0.500000000000
    (2578,   2,  -96):     0,   # 0xffffffa0,  frac=0.499981103553
    (2663,  82,  289):     1,   # 0x00000121,  frac=0.499973992497
    (2797, 110,  209):     0,   # 0x000000d1,  frac=0.500096974966
    (3241,  35,   86):     1,   # 0x00000056,  frac=0.500000000000
    (3263,   2,   76):     0,   # 0x0000004c,  frac=0.500000000000
    (3263,   9,  -14):     0,   # 0xfffffff2,  frac=0.500000000000
    (3263,  17,   -4):     0,   # 0xfffffffc,  frac=0.500000000000
    (3263,  18,  -24):     0,   # 0xffffffe8,  frac=0.500000000000
    (3263,  27,  -90):     1,   # 0xffffffa6,  frac=0.500000000000
    (3263,  30,   26):     0,   # 0x0000001a,  frac=0.500000000000
    (3298,   0,   95):    -1,   # 0x0000005f,  frac=0.500000000000
    (3316,   0,  110):    -1,   # 0x0000006e,  frac=0.500000000000
    (3502,  30,  127):     0,   # 0x0000007f,  frac=0.500000000000
    (3602,   0,   28):    -1,   # 0x0000001c,  frac=0.500000000000
    (3705,   0, -122):     0,   # 0xffffff86,  frac=0.500000000000
    (3876,   5,  -43):     0,   # 0xffffffd5,  frac=0.500000000000
    (3876,  10,   -8):     0,   # 0xfffffff8,  frac=0.500000000000
    (3876,  12,  -53):     0,   # 0xffffffcb,  frac=0.500000000000

    # For output = work/ANN.dTs.GHCN.CL.3
    ( 323,   0,  -17):    -1,   # 0xffffffef,  frac=0.500000000000
    ( 625,   0,   71):     0,   # 0x00000047,  frac=0.500000000000
    ( 908,   3,  -33):     0,   # 0xffffffdf,  frac=0.500000000000
    ( 908,   4,  -63):     0,   # 0xffffffc1,  frac=0.500000000000
    ( 908,   9,    2):     0,   # 0x00000002,  frac=0.500000000000
    ( 908,  11,  -63):     0,   # 0xffffffc1,  frac=0.500000000000
    (1006,   0,  -19):     0,   # 0xffffffed,  frac=0.500000000000
    (1007,   2,    8):     0,   # 0x00000008,  frac=0.500000000000
    (1007,   7,   33):     0,   # 0x00000021,  frac=0.500000000000
    (1007,  11,  -17):    -1,   # 0xffffffef,  frac=0.500000000000

    # For output = work/ANN.dTs.GHCN.CL.5
    (  96,  25,  -47):     0,   # 0xffffffd1,  frac=0.499925909420
    ( 232,   1,  173):     0,   # 0x000000ad,  frac=0.500017127746
    ( 339,   0,  107):     0,   # 0x0000006b,  frac=0.499802146304
    ( 360,  66,  -58):     0,   # 0xffffffc6,  frac=0.500240385412
    ( 408,   1,   91):    -1,   # 0x0000005b,  frac=0.500018695593
    ( 408,   3,   31):    -1,   # 0x0000001f,  frac=0.500018695593
    ( 408,  11,   21):    -1,   # 0x00000015,  frac=0.500018695593
    ( 408,  15,   56):    -1,   # 0x00000038,  frac=0.500018695593
    ( 408,  18,  -89):    -1,   # 0xffffffa7,  frac=0.499981304407
    ( 408,  28,   16):    -1,   # 0x00000010,  frac=0.500018695593
}


def fixRoundingIssue(fixData, yIdx, v, frac):
    recordIdx, fname = fixData
    key = recordIdx, yIdx, v
    hv = v
    if v < 0:
        hv = 0x100000000 + v
    adjust = _fixMap.get(key, 0)
    if key not in _fixMap:
        print "Potential round error:    (%4d, %3d, %4d): %5d,   # 0x%08x,  frac=%.12f" % (
                recordIdx, yIdx, v, adjust, hv, frac)
    return v + _fixMap.get(key, 0)


def annav(mon, nyrs, iy1, ibad, m1, emuBug=False, fixData=()):
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

            # Fortran appears to use single precision floats, but python uses
            # double precision. This can cause a fair number of differences
            # following summing and rounding. The following maths is to identify
            # potential value that suffer in this way, so they can (for now) be
            # adjusted to match the fortran results.
            absV = abs(v)
            w, frac = divmod(absV, 1)
            w = w + 0.5
            a, b = w * 0.9999950, w * 1.0000050
            if a < absV < b:
                iann[yIdx] = fixRoundingIssue(fixData, yIdx, iann[yIdx], frac)

            iy2n = iy1 + yIdx
            if iy1n == ibad:
                iy1n = iy2n

    return iy1n, iy2n, iann
    

def do_th_stuff():
    fname = "work/fort.2"
    f2 = fort.open("work/fort.2", "rb")
    f3 = ccc_binary.BufferedOutputRecordFile("work/fort.3")
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
            50, fmt="File %s: %%d records processed\n" % "fort.2")
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

        i1, i2, iann = annav(inRec.idata[0:nMonths], nyrs, iy1, ibad, m1,
                emuBug=(m2 % 12 == 0), fixData=(recordCount, fname))
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

        if options.verbose:
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
    do_th_stuff()


if __name__ == "__main__":
    import optparse
    usage = "usage: %prog [options]"
    parser = script_support.makeParser(usage)
    options, args = script_support.parseArgs(parser, __doc__, (0, 0))
    # Psyco simply does not play nicely with this program.
    # if not options.no_psyco:
    #     script_support.enablePysco(__file__, main,
    #            yearMonthEnumerate, fixRoundingIssue, annav, do_th_stuff)
    main(args)


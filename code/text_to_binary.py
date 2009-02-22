#!/usr/bin/env python
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
__docformat__ = "restructuredtext"

import sys
import array
from itertools import repeat

import fort
import ccc_binary
import script_support


def main(args):
    """The main for this module.

    This (currently) closely follows the form of the original fortran code in
    ``code/STEP2/text_to_binary.f``.
    """
    iyear1, lim, multm, ibad = 1880, 20, 2, 9999
    monthly, iok = [0]*12, [0]*12
    lat, lon, id, ht, hto, mmax = 0, 0, [0]*(multm + 1), 0, 0, [0]*(multm + 1)
    name, nameo, line, li, sid, sidi = "", "", "", "", "", ""


    lastyr = script_support.parseIntArg(0, args)
    MTOT = 12 * (lastyr - iyear1 + 1)
    if options.verbose >= 1:
        print "Last year with data: %s" % (lastyr,)
    idata = [array.array("i", [0] * MTOT) for i in range(multm + 1)]

    header = ccc_binary.CCHeader()
    header.title = 'GHCN V2 Temperatures (.1 C)'.ljust(80)
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
        if ict % 1000 == 0 and options.verbose >= 2:
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
                        if options.verbose >= 1:
                            f88.write("Dropped %s %s\n" % (id[m], nameo))

            mult = 0
            id1o = id1

        id[mult] = idfull
        idata[mult] = array.array("i", repeat(ibad, MTOT))
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
                    monmin = min(monmin, ix + m)
                    monmax = ix + m

            # F: goto 20
        # F: 30    continue

        f99.write("%s %s %s %s %s %s\n" % (
                mult, id[mult], mmax[mult],lim, monmin,monmax))
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
                if options.verbose >= 1:
                    f88.write("Dropped %s %s\n" % (id[m], nameo))

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


if __name__ == "__main__":
    usage = "usage: %prog [options] lastyear"
    parser = script_support.makeParser(usage)
    options, args = script_support.parseArgs(parser, __doc__, (1, 1))
    if not options.no_psyco:
        script_support.enablePysco(__file__, main, writeRecord)
    main(args)

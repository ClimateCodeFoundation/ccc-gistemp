#!/usr/bin/env python
"""Python replacement for code/STEP2/flags.f

Input files:

    work/fort.78

Output files:
    work/fort.2

The input file is a text file, with a single record per line. This is created
using the ``PApars.py`` script. This script adds the flags value at the end of 
each record, based on the record's contents. (Later the output file is given
a header line to create ``work/PApars.list``.)

"""
__docformat__ = "restructuredtext"


import ccc_binary
import script_support


# Constants
lshort = 7
slplim = 1.0
slpx = 0.5


def readRec(f):
    """Read a single record from the input file.

    The input file is text, line oriented, one record per line.

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
    return cc, IDc, sl1, sl2, knee, Yknee, sl0, Ymid,RMS,RMSl, iy1, iy2, iy1e, iy2e, l


def main(args):
    """The main for this module.

    """
    f78 = open("work/PApars.pre-flags")
    f2 = open("work/PApars.list", "w")
    f2.write('CCdStationID  slope-l  slope-r knee    Yknee    slope     Ymid      RMS     RMSl 3-rur+urb ext.range flag\n')
    nstat = sumch = nyrs = nstap = nyrsp = sumchp = 0
    nsl1 = nsl2 = ndsl = nswitch = nsw0 = nok = nokx = nsta = 0
    nshort = 0
    while True:
        data = readRec(f78)
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
        if knee < iy1+lshort or knee > iy2 - lshort:
            iflag += 1
        if knee < iy1+lshort or knee > iy2 - lshort:
            nshort += 1
        if abs(sl1) > slplim:
            iflag += 20
        if abs(sl1) > slplim:
            nsl1 += 1
        if abs(sl2) > slplim:
            iflag += 10
        if abs(sl2) > slplim:
            nsl2 += 1
        if abs(sl2 - sl1) > slplim:
            iflag += 100
        if abs(sl2 - sl1) > slplim:
            ndsl += 1
        # TODO: The small offset fixes a Fortran/Python rounding issue, but I
        #       am not happy with it. PAO, Mon Feb 23 2009.
        if abs(sl2 - sl1) + 0.00000001 > slpx:
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

    print " %-10s %4d %10.7f     %10.7f" % ("all", nsta,-sumch/nsta,-10.*sumch/nyrs)
    print " %-11s %8d  %10.7f    %10.7f" % ("urb warm", nstap,-sumchp/nstap,-10.*sumchp/nyrsp)
    print " %-11s %11d %11d %11d %11d %11d %11d" % (
        "# short,sl1,sl2,dsl,ok", nshort,nsl1,nsl2,ndsl,nok,nokx)
    print " %-11s  %11d %11d" % ("switches: all , else ok", nswitch,nsw0)


if __name__ == "__main__":
    import optparse
    usage = "usage: %prog [options]"
    parser = script_support.makeParser(usage)
    options, args = script_support.parseArgs(parser, __doc__, (0, 0))
    main(args)


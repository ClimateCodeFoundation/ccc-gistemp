#!/usr/bin/env python
"""Python replacement for code/STEP2/padjust.f

Input files:

    work/fort.1

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

def verbose(level, s):
    if options and level <= options.verbose:
        print >>sys.stderr, s


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


def main(args):
    """The main for this module.

    """
    f1 = open("work/PApars.list")
    header = f1.readline()
    # Read the first adjustment record.
    cc, IDc, sl1, sl2, knee, sl0, iy1, iy2, iy1e, iy2e, flag = readRec(f1)

    for iin in range(1, 7):
        inF = fort.open("work/Ts.GHCN.CL.%d" % iin)
        outF = fort.open("work/Ts.GHCN.CL.PA.%d" % iin, "wb")

        h = ccc_binary.CCHeader(inF.readline())
        ibad = h.info[6]
        if h.info[0] == ibad:
            continue
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
                    print " station   %9s  %s skipped" % (
                            rec_o.ID, rec_o.name)
                    continue

            else:
                # Record needs adjusting. Do so and get the next adjustment
                # entry.
                a, b = adj(h, rec_o.idata, sl1, sl2, knee, sl0, iy1e, iy2e, iy1, iy2,
                        flag, m1o, m2o)
                print " station   %9s  %s adjusted %s %s" % (
                        rec_o.ID, rec_o.name, (m1o, m2o), (a, b))
                aa = a - m1o
                bb = b - a + 1
                print "ADJ-1", rec_o.name, rec_o.count, len(rec_o.idata), aa, bb,
                rec_o.idata[:] = rec_o.idata[aa:aa + bb]
                rec_o.count = bb
                print rec_o.count, len(rec_o.idata)
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
                print " station   %9s  %s saved %11d" % (
                        rec_o.ID, rec_o.name, iin)
                break

            rec = ccc_binary.CCRecord(rec_o.m2 - rec_o.m1 + 1, data=s)
            if rec.ID != IDc:
                if rec.name[30:32] == " R" or rec.name[30] == "1":
                    outF.writeline(rec_o.binary)
                    print " station   %9s  %s saved %11d" % (
                            rec_o.ID, rec_o.name, iin)
                else:
                    rec_o.m1, rec_o.m2 = rec.m1, rec.m2
                    print " station   %9s  %s skipped" % (
                            rec.ID, rec.name)
                    continue

            else:
                a, b = adj(h, rec.idata, sl1, sl2, knee, sl0, iy1e, iy2e, iy1, iy2,
                        flag, rec_o.m1, rec_o.m2)
                print " station   %9s  %s adjusted %s %s" % (
                        rec.ID, rec.name, (rec_o.m1, rec_o.m2), (a, b))
                aa = a - rec_o.m1
                bb = b - a + 1
                print "ADJ-2", rec.name, rec.count, len(rec.idata), aa, bb,
                rec.idata[:] = rec.idata[aa:aa + bb]
                rec.count = bb
                print rec.count, len(rec.idata)
                rec_o.m1, rec_o.m2 = a, b
                outF.writeline(rec_o.binary)
                print " station   %9s  %s saved %11d" % (
                        rec_o.ID, rec_o.name, iin)
                IDc = -9999
                ddd = readRec(f1)
                if ddd[0] is not None:
                    cc, IDc, sl1, sl2, knee, sl0, iy1, iy2, iy1e, iy2e, flag = ddd

            rec_o = rec.copy()


import sys
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


if __name__ == "__main__":
    import optparse
    usage = "usage: %prog [options]"
    parser = script_support.makeParser(usage)
    options, args = script_support.parseArgs(parser, __doc__, (0, 0))
    main(args)


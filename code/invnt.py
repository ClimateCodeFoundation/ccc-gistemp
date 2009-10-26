#!/usr/bin/env python
"""Python replacement for code/STEP2/invnt.f

This basically scans the files ``Ts.GHCN.CL.1``, ``Ts.GHCN.CL.2``, etc.
(created by ``trim_binary``) and extracts station information.

As used in practice, the output is written to ``Ts.GHCN.CL.station.list``.
"""
__docformat__ = "restructuredtext"


# Implementation notes:
#
# The fortran program contains a fair amount of unused code (probably used
# for debugging). None of that is replicated here.

import sys
import os

import fort
import ccc_binary
import script_support


iyrbeg = 1880
iyrend = 3001
ID = 0
nrecpy = 12
ndpryr = nrecpy
monm = ndpryr * (iyrend - iyrbeg + 1)

def invntSingleFile(fileNumber):
    fname = "%s.%d" % (args[0], fileNumber)
    progress = script_support.countReport(
            50, fmt="File %s: %%d records processed\n" % fname)
    f15 = fort.open(fname, "rb")
    header = ccc_binary.CCHeader(data=f15.readline())
    n1 = header.info[0]
    n2 = header.info[8]

    if header.info[6 - 1] < iyrbeg:
        sys.exit("IYRBEG TOO HIGH")
    if header.info[4 - 1] > monm:
        sys.exit("IYREND TOO LOW %s > %s" % (header.info[4 - 1], monm))
    mbad = header.info[7 - 1]

    if header.info[1 - 1] == mbad:
        # There is no data
        return False

    for recordCount, s in enumerate(f15):
        length = n2 - n1 + 1
        sys.stdout.flush()
        rec = ccc_binary.CCRecord(length, data=s)
        if ID <= rec.ID:
            sys.stdout.write(
                "%s %9d %-30s lat,lon (.1deg)%5d%6d %s%s%s cc=%3s\n" % (
                    fname, rec.ID, rec.name[:30], rec.Lat, rec.Lon,
                    rec.name[31], rec.name[30], rec.name[32],
                    rec.name[33:36],
                ))
            n1 = rec.m1
            n2 = rec.m2

            if options.verbose:
                progress.next()

    return True


def main(args):
    """The main for this module.

    This simply invokes invntSingleFile for n = 1, 2, ... 6.
    """
    for n in range(1, 7):
        if not invntSingleFile(n):
            # No data?
            break


if __name__ == "__main__":
    import optparse
    usage = "usage: %prog [options]"
    parser = script_support.makeParser(usage)
    options, args = script_support.parseArgs(parser, __doc__, (1, 1))
    main(args)


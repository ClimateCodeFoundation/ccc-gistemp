#!/usr/bin/env python
"""Python replacement for code/STEP2/trim_binary.f

The program reads the input files ``work/Ts.bin1``, ``work/Ts.bin2``, etc. and
trims leading an trailing blocks of month data, where ``abs(md) > 8000``. The
reduced size records are writtent o mathcing files ``work/Ts.bin1.trim``,
``work/Ts.bin2.trim``, etc.
"""
__docformat__ = "restructuredtext"

import sys

import fort
import ccc_binary
import script_support


def trimSingleFile(fileNumber):
    """The main for this module.

    Note:
        This (currently) tries to closely follow the form of the original
        fortran code in ``code/STEP2/trim_binary.f``. Hence it is too long, too
        nested, etc compared to typical Python code.

    """
    filei = "work/Ts.bin%s" % fileNumber
    f2 = fort.open(filei, "rb")
    fileo = "work/Ts.GHCN.CL.%s" % fileNumber
    f3 = ccc_binary.BufferedOutputRecordFile(fileo)

    header = ccc_binary.CCHeader(f2.readline())
    f3.writeRecord(header)
    if options.verbose >= 3:
        print "Input file %r header details:" % filei
        print "    %s" % header.title.rstrip()
        for i, v in enumerate(header.info):
            print "    Info[%d] = %10d / %08x" % (i, v, v)

    i4 = header.info[3]
    inRec = ccc_binary.CCRecord(i4)
    marker = header.info[6]

    # TODO: What if m1 == m2 == marker. Should that be:
    #
    # - a fatal condition.
    # - cause an empty record o be written.
    # - cause no record to be written.

    progress = script_support.countReport(
            50, fmt="File %d: %%d records processed\n" % fileNumber, f=sys.stdout)
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

        if options.verbose:
            progress.next()


def main(args):
    """The main for this module.

    This simply invokes trimSingleFile for n = 1, 2, ... 6.
    """
    for n in range(1, 7):
        trimSingleFile(n)


if __name__ == "__main__":
    import optparse
    usage = "usage: %prog [options]"
    parser = script_support.makeParser(usage)
    options, args = script_support.parseArgs(parser, __doc__, (0, 0))
    main(args)


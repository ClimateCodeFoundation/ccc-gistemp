#!/usr/bin/env python
"""Python replacement for code/STEP2/split_binary.f

The program reads the input file ``work/Ts.bin`` and splits it to produce
6 files ``Ts.bin1``, ``Ts.bin2``, ... ``Ts.bin6``.
"""
__docformat__ = "restructuredtext"

import sys

import fort
import ccc_binary
import script_support


def main(args):
    """The main for this module.

    This (currently) tries to closely follow the form of the original fortran
    code in ``code/STEP2/split_binary.f``.
    """
    inPath = "work/Ts.bin"
    try:
        inF = fort.open(inPath)
    except IOError, exc:
        sys.exit("Open failed for %s\n%s" % (inPath, exc))
    header = ccc_binary.CCHeader(data=inF.readline())
    if options.verbose >= 1:
        print "Input file header details:"
        print "    %s" % header.title.rstrip()
        for i, v in enumerate(header.info):
            print "    Info[%d] = %10d / %08x" % (i, v, v)

    # Open each output file and write the header.
    outF = {}
    for n in range(51, 56+1):
        fileo = "work/Ts.bin%d" % (n - 50,)
        outF[n] = fort.open(fileo, "wb")
        outF[n].writeline(header.binary)

    # Process the input file, copying each records to one of the the 6 output
    # files.
    rec = ccc_binary.CCRecord(header.info[3])
    for recordCount, s in enumerate(inF):
        rec.setFromBinary(s)
        mu = rec.Lat + 899
        if mu < 0:
            mu = 0
        mu = 56 - mu // 300
        outF[mu].writeline(rec.binary)

        if options.verbose >= 2:
            print "Write record %d to %d" % (recordCount, mu - 50)

    if options.verbose >= 1:
        print "Processed records =", recordCount


if __name__ == "__main__":
    import optparse
    usage = "usage: %prog [options]"
    parser = script_support.makeParser(usage)
    options, args = script_support.parseArgs(parser, __doc__, (0, 0))
    main(args)

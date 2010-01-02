#!/usr/bin/env python
"""Python replacement for code/STEP4_5/trimSBBX.f

Input files:

    tbd

Output files:

    tbd

"""
__docformat__ = "restructuredtext"


import sys
import struct

import fort
import ccc_binary
import script_support


km0 = 12
monm0 = km0 * (3000 - 1850 + 1)
nsbbx = 8000


class Rec(ccc_binary.CCBin):
    """Holds input file records.

    """
    def __init__(self, nAvg, nLat, **kwargs):
        self.nAvg = nAvg
        self.nLat = nLat
        super(Rec, self).__init__(**kwargs)

    @property
    def _format(self):
        return "%df%di" % (self.nAvg, self.nLat)

    def initEmpty(self):
        self.avg = [0] * max(self.nAvg, monm0)
        self.lat = [0] * max(self.nLat, 7)

    def getElements(self):
        return self.avg[:self.nAvg] + self.lat[:self.nLat]

    def setFromElements(self, elements):
        self.initEmpty()
        self.avg[:self.nAvg] = list(elements[:self.nAvg])
        self.lat[:self.nLat] = list(elements[self.nAvg:self.nAvg + self.nLat])


class TrimmedRec(Rec):
    def __init__(self, nAvg, nLat, **kwargs):
        super(TrimmedRec, self).__init__(nAvg, nLat, **kwargs)
        self.mln = 0

    @property
    def _format(self):
        return "i7i%df" % (self.nAvg, )

    def getElements(self):
        return [self.mln] + self.lat[:7] + self.avg[:self.nAvg]

    def set_from_rec(self, rec):
        self.nAvg = rec.nAvg
        self.nLat = rec.nLat
        self.avg = list(rec.avg)
        self.lat = list(rec.lat)


def findl6(avg, ml, xbad):
    """Find the value for element 6 of the LAT array.

    This counts the number of values in ``avg[:ml]`` that are good (i.e.
    not equal to ``xbad``). This is the return value and is the value to
    be used for the 6th element in the ``LAT`` array in an output record.

    :Return:
        The number of good average values.

    """
    return len([el for el in avg[:ml] if el != xbad])


def main(in_path, out_path):
    """The main for this module.

    """
    f10 = fort.open(in_path)
    f10.bos = ">"
    f11 = fort.open(out_path, "wb")
    f11.bos = ">"
    h = ccc_binary.SBBX_Header(f10.readline(), bos=">")

    km = 1
    if h.mavg == 6:
        km = 12
    if h.mavg == 7:
        km = 4
    # Here Fortran code verifies that the input file gives km=12, but that will
    # always be true for the CCC code.
    if km != km0:
        sys.stderr.write("ERROR: CHANGE KM0\n")
        return 1

    ml = h.nm
    if monm0 < ml:
        sys.stderr.write("Set monm0 at least to %s\n" % ml)
        return 1

    mbad = last = h.bad
    xbad = float(mbad)
    newHdr = h.copy()
    newHdr.set_bos(">")
    newHdr.recsize = h.nm + 8
    nLat = h.recsize - h.nm

    rec = Rec(ml, nLat, bos=">")
    reco = TrimmedRec(ml, nLat, bos=">")
    rec.lat[4:6] = [0, 1, 0]
    reco.lat[4:6] = [0, 1, 0]
    temp_rec = Rec(ml, nLat, bos=">")
    temp_rec.setFromBinary(data=f10.readline())
    reco.set_from_rec(temp_rec)
    newHdr.ml = 1

    if nLat == 4:
        reco.lat[5] = findl6(reco.avg, ml, xbad)
    # If first SBBX has data then set ``newHdr.ml`` to ``ml``.
    if reco.lat[5] > 0:
        newHdr.ml = ml

    f11.writeline(newHdr.binary)

    for n in range(2, nsbbx + 1):
        rec = Rec(ml, nLat, data=f10.readline(), bos=">")
        reco.mln = 1
        reco.nAvg = ml
        if nLat == 4:
            rec.lat[5] = findl6(rec.avg, ml, xbad)
        if rec.lat[5] > 0:
            reco.mln = ml
        if reco.lat[5] == 0:
            data = tuple([reco.mln] + reco.lat + [xbad])
            bin = struct.pack(">i7if", *data)
            f11.writeline(bin)
        else:
            f11.writeline(reco.binary)
        reco.set_from_rec(rec)
        reco.nAvg = ml

    reco.mln = 0
    if reco.lat[5] == 0:
        data = tuple([reco.mln] + reco.lat + [xbad])
        bin = struct.pack(">i7if", *data)
        f11.writeline(bin)
    else:
        f11.writeline(reco.binary)
    return 0


if __name__ == "__main__":
    import optparse
    usage = "usage: %prog [options] [infile] [outfile]"
    usage += "\n\nThe infile, outfile default to work/fort.10, work/fort.11"
    parser = script_support.makeParser(usage)
    options, args = script_support.parseArgs(parser, __doc__, (0, 2))
    sys.exit(main(*args))

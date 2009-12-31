#!/usr/bin/env python
"""Python version of convert1.HadR2_mod4.py.

"""
__docformat__ = "restructuredtext"


import sys
import string
import struct
import array
import copy

import fort
import ccc_binary
import script_support


# Constants (from original fortran parameters)
IM = 360
JM = 180
IOFF = IM // 2
IYRBEG = 1880
IYREO = 3000
NYMAX = 15
MONX = 12 * (IYREO - IYRBEG + 1)


def make_3d_array(a, b, c):
    """Create an array with three dimensions.

    Actually a list-of-lists-of-lists, but the result can be treated as an
    array.

    """
    x = [0.0] * c
    arr = []
    for i in range(a):
        arr.append([list(x) for j in range(b)])

    return arr


def sread(f, na):
    data = f.readline()
    mnext, lts, ltn, lnw, lne, x, x, x = struct.unpack(">lllllfff", data[:32])
    a = struct.unpack(">%df" % na, data[32:])
    return mnext, lts, ltn, lnw, lne, a


def load_oisstv2_mod4_clim():
    f14 = fort.open("input/oisstv2_mod4.clim")
    try:
        f14.bos = ">"
        data = f14.readline()
    finally:
        f14.close()
    clim_title, data = data[:80], data[80:]

    clim = make_3d_array(IM, JM, 12)
    ii = jj = kk = 0
    for p in range(0, len(data), 4):
        v, = struct.unpack(">f", data[p:p+4])
        int_v, = struct.unpack("<L", data[p:p+4])
        clim[ii][jj][kk] = v
        ii += 1
        if ii >= IM:
            ii = 0
            jj += 1
            if jj >= JM:
                jj = 0
                kk += 1

    return clim


def load_sst_data(iyr1, mo1, mo2):
    # Read in the SST data for recent years (unit 11)
    sst = make_3d_array(IM, JM, 12 * NYMAX)

    for iyr in range(iyr1, iyr1 + 1):
        for mon in range(mo1, mo2 + 1):
            path = "input/oiv2mon.%04d%02d" % (iyr, mon)
            sys.stdout.write(" trying to read %-70s\n" % path)
            f11 = fort.open(path)
            f11.bos = ">"
            data = f11.readline()
            iyr0, imo = struct.unpack(">LL", data[:8])
            iyre = iyr
            moe = mon
            if iyr != iyr0:
                sys.stderr.write("years not ok %s %s\n" % (iyr, iyr0))
            if mon != imo:
                sys.stderr.write("months not ok %s %s\n" % (mon, imo))

            data = f11.readline()
            f11.close()
            kk = 12 * (iyr - iyr1) + mon - 1
            ii = jj = 0
            for p in range(0, len(data), 4):
                v, = struct.unpack(">f", data[p:p+4])
                int_v, = struct.unpack("<L", data[p:p+4])
                sst[ii][jj][kk] = v
                ii += 1
                if ii >= IM:
                    ii = 0
                    jj += 1

    return sst, moe, iyre


def open_curr_SBBX():
    f1 = fort.open("input/SBBX.HadR2")
    f1.bos = ">"
    inH = ccc_binary.SBBX_Header(f1.readline(), bos=f1.bos)
    print inH.title

    return f1, inH


def main(args):
    f1, inH = open_curr_SBBX()

    bad = inH.bad
    xbad = float(bad)
    mnow = inH.m1

    if inH.yrbeg != IYRBEG:
        sys.stderr.write("Error: iyrbeg-new %s iyrbeg-old %s\n" % (
            IYRBEG, inH.yrbeg))
        sys.exit("Error: iyrbeg inconsistent")
    sys.stdout.write(" SBBX.HadR2  opened\n")

    # Paul: I do not know what this file is for. The fortran code opens the
    # file but does not seem to use it.
    # f24 = fort.open("input/SBBX_LtSN.LnWE.dat")
    # f24.bos = ">"

    iyr1, mo1, mo2 = [int(v) for v in args]
    clim = load_oisstv2_mod4_clim()
    sst, moe, iyre = load_sst_data(iyr1, mo1, mo2)

    outH = ccc_binary.SBBX_Header(bos="<")
    for attr in "m1 kq mavg nm recsize yrbeg bad trace".split():
        setattr(outH, attr, getattr(inH, attr))
    outH.m1 = 1                        # output file NOT trimmed/reorganized
    outH.nm = 12 * (iyr1 - IYRBEG + 1) # wipe out data after new data
    outH.recsize = outH.nm + 4         # old non-trimmed format
    outH.title = inH.title[:40] + " Had: 1880-11/1981, oi2: 12/1981-%2d/%04d" % (
            moe, iyre)

    f2 = fort.open("work/SBBX.HadR2.upd_py", "wb")
    f2.writeline(outH.binary)

    #
    # Interpolate to Sergei's subbox grid
    #
    to_ = [0.0] * MONX
    moff = (iyr1 - IYRBEG) * 12
    for n in range(8000):
        for m in range(MONX):
            to_[m] = bad
        mnext, lts, ltn, lnw, lne, a = sread(f1, mnow)
        to_[:mnow] = a
        mnow = mnext
        js = (18001 + (lts + 9000) * JM)/18000
        jn = (17999 + (ltn + 9000) * JM)/18000
        iw = (36001 + (lnw + 18000) * IM)/36000 + IOFF
        ie = (35999 + (lne + 18000) * IM)/36000 + IOFF
        if ie > IM:
            iw = iw - IM
            ie = ie - IM
        if iw > IM:
            sys.exit('Fatal: iw>im')
        if ie > IM:
            sys.exit('Fatal: ie>im')
        if iw < 1:
            sys.exit('Fatal: iw<1')
        if ie < 1:
            sys.exit('Fatal: ie<1')

        for m in range(mo1, mo2 + 1):
            month = (m - 1) % 12
            kt = 0
            ts = 0.0
            for j in range(js -1, jn):
                for i in range(iw - 1, ie):
                    if sst[i][j][m-1] < -1.77 or clim[i][j][month] == xbad:
                        continue
                    kt += 1
                    ts += sst[i][j][m-1] - clim[i][j][month]
            
            to_[m + moff - 1] = bad
            if kt > 0:
                to_[m + moff - 1] = ts / kt

        b = outH.nm
        data = to_[:b]
        data.extend([lts, ltn, lnw, lne])
        s = struct.pack("%dfllll" % (b,), *data)
        f2.writeline(s)


if __name__ == "__main__":
    import optparse
    usage = "usage: %prog [options] year month1 month2"
    parser = script_support.makeParser(usage)
    options, args = script_support.parseArgs(parser, __doc__, (3, 3))
    main(args)

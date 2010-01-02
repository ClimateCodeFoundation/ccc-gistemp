#!/usr/bin/env python
"""Python version of convert1.HadR2_mod4.py.

Without arguments this will look for all the ``oiv2mon.*`` files in the
``input_files`` directory and use them to buid a new ``input/SBBX.HadR2`` file.
Alternatively, you can provide 2 or 3 arguments to select the set of input
files. The arguments can be in one of two forms::

    year month1 month2
    year1 year2 [month2]

The first form can be used to select a number of months in a given year. For
example, ``2009 1 5``, will select Jan-May in 2009.

The second form is typically used to select a number of complete years. For
example, ``2005 2008``, will select all the months from Jan 2005 to Dec 2008.
If the optional ``month2`` is provided, then only part of the final year is
selected. For example, ``2005 2009 5``, will select Jan 2005 to May 2009.

"""
__docformat__ = "restructuredtext"


import os
import sys
import struct
import re
from glob import glob
try:
    import gzip
except ImportError:
    # No gzip module, this is not critical at this time.
    gzip = None

import fort
import ccc_binary
import script_support
from script_support import trace
import trimSBBX


# This supports using this script as a module.
class _Struct(object):
    pass
options = _Struct()
options.verbose = 0


# Constants (from original fortran parameters)
IM = 360
JM = 180
IOFF = IM // 2
IYRBEG = 1880
IYREO = 3000
NYMAX = 15
MONX = 12 * (IYREO - IYRBEG + 1)


# This is used to extract the end month/year from the title of the SBBX file.
rTitle = re.compile(r"Monthly Sea Surface Temperature anom \(C\)"
        " Had: 1880-11/1981, oi2: 12/1981- *(?P<month>\d+)/(?P<year>\d+)")


def make_3d_array(a, b, c):
    """Create an array with three dimensions.

    Actually a list-of-lists-of-lists, but the result can be treated as an
    array with dimentions ``[a][b][c]``.

    """
    x = [0.0] * c
    arr = []
    for i in range(a):
        arr.append([list(x) for j in range(b)])

    return arr


def sread(f, na):
    data = f.readline()
    mnext, lts, ltn, lnw, lne, x, x, x = struct.unpack(">lllllfff", data[:32])
    try:
        a = struct.unpack(">%df" % na, data[32:])
    except MemoryError:
        print "!!!!", len(data[32:]), na, hex(na)
        raise
    except struct.error:
        print "!!!!", len(data[32:]), na
        raise
    return mnext, lts, ltn, lnw, lne, a


def load_oisstv2_mod4_clim():
    src = "input_files/oisstv2_mod4.clim"
    if not os.path.exists("input_files/oisstv2_mod4.clim"):
        sys.stderr.write("File %r is missing\n" % src)
        return None

    f14 = fort.open(src)
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
        #int_v, = struct.unpack(">L", data[p:p+4])
        clim[ii][jj][kk] = v
        ii += 1
        if ii >= IM:
            ii = 0
            jj += 1
            if jj >= JM:
                jj = 0
                kk += 1

    return clim


def load_sst_data(iyr1, n_years, dates):
    # Read in the SST data for recent years (unit 11)
    sst = make_3d_array(IM, JM, 12 * n_years)

    trace(0, " Reading input_files/oiv2mon... files.")
    for iyr, mon in dates:
        path = "input_files/oiv2mon.%04d%02d" % (iyr, mon)
        trace(1, " trying to read %-70s" % path)
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
            # int_v, = struct.unpack(">L", data[p:p+4])
            sst[ii][jj][kk] = v
            # if v != 9999.0:
            #     print "IN:", ii+1, jj+1, kk+1, v, hex(int_v)
            ii += 1
            if ii >= IM:
                ii = 0
                jj += 1

    return sst, moe, iyre


def open_curr_SBBX():
    f1 = fort.open("input/SBBX.HadR2")
    f1.bos = ">"
    inH = ccc_binary.SBBX_Header(f1.readline(), bos=f1.bos)
    if options and options.verbose >= 3:
        for attr in "ml kq mavg nm recsize yrbeg bad trace".split():
            trace(3, "  %-5s = %s", attr, getattr(inH, attr))

    return f1, inH


def decompress_gz(p):
    dp = p[:-3]
    if os.path.exists(dp):
        return 2
    if gzip is None:
        return 1
    z = gzip.GzipFile(p, "rb")
    f = open(dp, "wb")
    while True:
        block = z.read(4096)
        if not block:
            break
        f.write(block)
    f.close()
    z.close()
    return 0


def parse_args(args):
    try:
        args = [int(v) for v in args]
    except ValueError:
        sys.exit("Arguments must be integers")

    if len(args) == 2:
        y1, y2 = args
        for y in args:
            if not 1981 <= y <= IYREO:
                sys.exit("Years must be between 1981 and %s" % IYREO)
        y1, y2 = args
        if y2 < y1:
            sys.exit("Second year cannot be before first year")
        return y1, 1, y2, 12

    if len(args) == 3:
        a, b, c = args
        if 1 <= b <= 12:
            # Assume args of the form ``year month1 month2``
            y1 = a
            m1, m2 = b, c
            if not 1981 <= y1 <= IYREO:
                sys.exit("Years must be between 1981 and %s" % IYREO)
            for m in (m1, m2):
                if not 1 <= m <= 12:
                    sys.exit("Months must be between 1 and 12")
            if m2 < m1:
                sys.exit("Second month cannot be before first month")
            return y1, m1, y1, m2

        else:
            # Assume args of the form ``year1 year2 month2``
            for y in (a, b):
                if not 1981 <= y <= IYREO:
                    sys.exit("Years must be between 1981 and %s" % IYREO)
                if not 1 <= c <= 12:
                    sys.exit("Months must be between 1 and 12")
            y1, m1, y2, m2 = a, 1, b, c
            if y2 < y1:
                sys.exit("Second year cannot be before first year")
            return y1, m1, y2, m2

    assert False, "This should never happen!"


def select_input_files(args, curr_moe=None, curr_yre=None,
        use_all_oiv2mon=False):
    dates = []

    if not os.path.isdir("input_files"):
        sys.stderr.write(" There is no input_files directory\n")
        sys.stderr.write(" The optional step 4 requires this directory\n")
        return 0, None, None, None

    # Try to ensure the input files are decompressed.
    gzipped_paths = glob(
            "input_files/oiv2mon.[0-9][0-9][0-9][0-9][0-9][0-9].gz")
    for i, p in enumerate(gzipped_paths):
        ret = decompress_gz(p)
        if ret == 1:
            sys.stderr.write("Warning: Unable to uncompress %s\n" % p)
            sys.stderr.write("         Please add python gzip support or"
                " unzip the oiv2mon files manually\n")
            return 0, None, None, None
        elif ret == 0:
            trace(1, " Ungzipped %s", p)

    if not args and use_all_oiv2mon:
        potential_paths = glob(
                "input_files/oiv2mon.[0-9][0-9][0-9][0-9][0-9][0-9]")
        for p in sorted(potential_paths):
            y, m = [int(s) for s in (p[-6:-2], p[-2:])]
            if not 1981 <= y <= IYREO or not 1 <= m <= 12:
                sys.stderr.write("Ignoring %r\n" % p)
                continue
            dates.append((y, m))

    elif not args and not use_all_oiv2mon:
        # Only select the files necessary to add data.
        start_year = curr_yre
        if curr_moe == 12:
            start_year += 1
        potential_paths = glob(
                "input_files/oiv2mon.[0-9][0-9][0-9][0-9][0-9][0-9]")
        for p in sorted(potential_paths):
            y, m = [int(s) for s in (p[-6:-2], p[-2:])]
            if y < start_year:
                continue
            if not 1981 <= y <= IYREO or not 1 <= m <= 12:
                sys.stderr.write("Ignoring %r\n" % p)
                continue
            dates.append((y, m))

        if not dates or dates[-1] <= (curr_yre, curr_moe):
            # There are no newer data.
            return None, None, None, None

    else:
        iyr1, mo1, iyr2, mo2 = parse_args(args)
        y, m = iyr1, mo1
        while (y, m) <= (iyr2, mo2):
            p = "input_files/oiv2mon.%04d%02d" % (y, m)
            if not os.path.exists(p):
                sys.stderr.write(
                        "Input file %r does not exist - ignoring\n" % p)
            else:
                dates.append((y, m))
            m += 1
            if m > 12:
                m = 1
                y += 1

    if not dates:
        sys.stderr.write(
                "No oiv2mon... files found\n")
        return 0, None, None, None

    dates.sort()
    n_years = dates[-1][0] - dates[0][0] + 1
    return dates[0][0], dates[-1][0], n_years, dates


def main(args=(), use_all_oiv2mon=False, verbose=None):
    if verbose is not None:
        options.verbose = verbose
        script_support.options.verbose = verbose

    f1, inH = open_curr_SBBX()

    bad = inH.bad
    xbad = float(bad)
    mnow = inH.ml

    if inH.yrbeg != IYRBEG:
        sys.stderr.write("Error: iyrbeg-new %s iyrbeg-old %s\n" % (
            IYRBEG, inH.yrbeg))
        sys.exit("Error: iyrbeg inconsistent")
    trace(1, " SBBX.HadR2 opened")
    trace(1, " Orig header: %s", inH.title)

    m = rTitle.match(inH.title)
    if m is None:
        sys.stderr.write("The title in %s does not look right\n"
                % "SBBX.HadR2")
        sys.stderr.write("Unable to determine end month/year from:\n")
        sys.stderr.write("  %r\n" % inH.title)
        sys.exit(1)
    curr_moe, curr_yre = [int(v) for v in m.groups()]

    iyr1, iyr2, n_years, dates = select_input_files(args, curr_moe, curr_yre,
            use_all_oiv2mon=use_all_oiv2mon)
    if iyr1 is None:
        sys.stdout.write(" The SBBX.HadR2 appears to be up to date,"
                         " leaving it unchanged\n")
        return 0
    elif iyr1 == 0:
        sys.stderr.write(" Unable to update SBBX.HadR2\n")
        return 0

    if 0:
        for y, m in dates:
            path = "input_files/oiv2mon.%04d%02d" % (y, m)
            print path

    sst, moe, iyre = load_sst_data(iyr1, n_years, dates)
    clim = load_oisstv2_mod4_clim()
    if clim is None:
        sys.stderr.write(" Unable to update SBBX.HadR2\n")
        return 1
    trace(1, " All input files read")

    outH = ccc_binary.SBBX_Header(bos=">")
    for attr in "ml kq mavg nm recsize yrbeg bad trace".split():
        setattr(outH, attr, getattr(inH, attr))
    outH.ml = 1                        # output file NOT trimmed/reorganized
    outH.nm = 12 * (iyr2 - IYRBEG + 1) # wipe out data after new data
    outH.recsize = outH.nm + 4         # old non-trimmed format
    outH.title = inH.title[:40] + " Had: 1880-11/1981, oi2: 12/1981-%2d/%04d" % (
            moe, iyre)

    f2 = fort.open("work/SBBX.HadR2", "wb")
    f2.bos = ">"
    f2.writeline(outH.binary)

    #
    # Interpolate to Sergei's subbox grid
    #
    trace(0, " Interpolating to Sergei's subbox grid")
    to_ = [0.0] * MONX
    moff = (iyr1 - IYRBEG) * 12
    for n in range(8000):
        if n and n % 500 == 0:
            trace(1, " %s of 8000 positions processed", n)
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

        for y, mm in dates:
            m = (y - iyr1) * 12 + mm
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
        s = struct.pack(">%dfllll" % (b,), *data)
        f2.writeline(s)

    # Close the updated file and trim.
    f2.close()
    trace(0, " Trimming the updated SBBX file")
    tmp_out = "work/SBBX.HadR2.trim"
    ret = trimSBBX.main("work/SBBX.HadR2", tmp_out)
    if ret == 0:
        # Trimming was successful. So rename the working file to become
        # the input/SBBX.HadR2
        trace(0, " Replacing input/SBBX.HadR2 with new file")
        os.rename(tmp_out, "input/SBBX.HadR2")


if __name__ == "__main__":
    import optparse
    usage =    "usage: %prog [options] [year month1 month2]"
    usage += "\n       %prog [options] [year1 year2 [month2]]"
    parser = script_support.makeParser(usage)
    parser.add_option("--use-all-oiv2mon", action="store_true",
        help="Use all the oiv2mon files in input_files")
    options, args = script_support.parseArgs(parser, __doc__, (0, 3))
    if len(args) not in (0, 2, 3):
        parser.error("Expected either zero, 2 or 3 arguments")
    sys.exit(main(args, options.use_all_oiv2mon))

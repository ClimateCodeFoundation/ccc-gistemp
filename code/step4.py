#!/usr/bin/env python
# $URL$
# $Rev$

"""Python version of convert1.HadR2_mod4.py.

Without arguments this will look for ``input/SBBX.HadR2`` and all the
``oiv2mon.*`` files in the ``input`` directory and use them to buid a
new ``work/SBBX.HadR2`` file.  Alternatively, you can provide 2 or 3
arguments to select the set of input files. The arguments can be in
one of two forms::

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


def load_oisstv2_mod4_clim():
    src = "input/oisstv2_mod4.clim"
    if not os.path.exists("input/oisstv2_mod4.clim"):
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

    for iyr, mon in dates:
        path = "input/oiv2mon.%04d%02d" % (iyr, mon)
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
            sst[ii][jj][kk] = v
            ii += 1
            if ii >= IM:
                ii = 0
                jj += 1

    return sst, moe, iyre


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


def select_input_files(curr_moe=None, curr_yre=None, args=[],
        use_all_oiv2mon=False):
    dates = []

    if not os.path.isdir("input"):
        sys.stderr.write(" There is no input directory\n")
        return 0, None, None, None

    # Try to ensure the input files are decompressed.
    gzipped_paths = glob(
            "input/oiv2mon.[0-9][0-9][0-9][0-9][0-9][0-9].gz")
    for i, p in enumerate(gzipped_paths):
        ret = decompress_gz(p)
        if ret == 1:
            sys.stderr.write("Warning: Unable to uncompress %s\n" % p)
            sys.stderr.write("         Please add python gzip support or"
                " unzip the oiv2mon files manually\n")
            return 0, None, None, None

    # Only select the files necessary to add data.
    start_year = curr_yre
    if curr_moe == 12:
        start_year += 1
    potential_paths = glob(
            "input/oiv2mon.[0-9][0-9][0-9][0-9][0-9][0-9]")
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

    if not dates:
        sys.stderr.write(
                "No oiv2mon... files found\n")
        return 0, None, None, None

    dates.sort()
    n_years = dates[-1][0] - dates[0][0] + 1
    return dates[0][0], dates[-1][0], n_years, dates


class Step4Iterator(object):
    def __init__(self, record_source):
        self.box_source = iter(self._check_for_new_data(record_source))
        self.meta = self.box_source.next()

    def _check_for_new_data(self, record_source):
        # The title string of the input source contains some metadata details.
        # So we parse the title to figure out which years and months the
        # source currently contains.
        m = rTitle.match(record_source.meta.title)
        if m is None:
            sys.stderr.write("The title in %s does not look right\n"
                    % "SBBX.HadR2")
            sys.stderr.write("Unable to determine end month/year from:\n")
            sys.stderr.write("  %r\n" % inH.title)
            sys.exit(1)
        curr_moe, curr_yre = [int(v) for v in m.groups()]

        # Look for 'oiv2mon...' files that contain additional data.
        iyr1, iyr2, n_years, dates = select_input_files(curr_moe, curr_yre)
        if iyr1 == 0:
            sys.stderr.write(" Unable to update SBBX.HadR2\n")
            return record_source

        elif iyr1 is None:
            sys.stdout.write(" The input SBBX.HadR2 appears to be up to date,"
                             " using it unchanged\n")
            return record_source

        # We have new data. Create an generator that provides the new data
        # set.
        return MergeReader(record_source, iyr1, iyr2, n_years, dates)

    def __iter__(self):
        return self._it()

    def _it(self):
        yield self.meta
        for box in self.box_source:
            yield box


def step4(data):
    """Step 4 of GISTEMP processing.  This is a little unusual compared
    to the other steps.  The effect of running this step is to update
    the SBBX.HadR2 file.  The input data is a pair of iterables.
    It produces a zipped output.
    
    """
    import itertools

    land,ocean = data
    oceanresult = Step4Iterator(ocean)
    return itertools.izip(land, oceanresult)


class MergeReader(object):
    def __init__(self, reader, iyr1, iyr2, n_years, dates):
        self.reader = iter(reader)
        self.iyr1 = iyr1
        self.dates = dates
        self.sst, moe, iyre = load_sst_data(self.iyr1, n_years, dates)
        self.clim = load_oisstv2_mod4_clim()
        if self.clim is None:
            sys.stderr.write(" Unable to update SBBX.HadR2\n")
            # TODO: What now?
            return
    
        self.meta = self.reader.next()
        self.meta.monm = 12 * (iyr2 - IYRBEG + 1)
        self.meta.monm4 = self.meta.monm + 8
        self.meta.title = (reader.meta.title[:40] 
                + " Had: 1880-11/1981, oi2: 12/1981-%2d/%04d" % (moe, iyre))

    def __iter__(self):
        return self.it()

    def it(self):
        yield self.meta

        #
        # Interpolate to Sergei's subbox grid
        #
        moff = (self.iyr1 - IYRBEG) * 12
        boxes = self.reader
        sst = self.sst
        clim = self.clim
        for box in boxes:
            box.pad_with_missing(self.meta.monm)
            lts = box.lat_S
            ltn = box.lat_N
            lnw = box.lon_W
            lne = box.lon_E

            js = (18001 + (lts + 9000) * JM)/18000
            jn = (17999 + (ltn + 9000) * JM)/18000
            iw = (36001 + (lnw + 18000) * IM)/36000 + IOFF
            ie = (35999 + (lne + 18000) * IM)/36000 + IOFF
            if ie > IM:
                iw = iw - IM
                ie = ie - IM
            if iw > IM:
                sys.exit('Fatal: iw>im, %s %s' % (iw,IM))
            if ie > IM:
                sys.exit('Fatal: ie>im')
            if iw < 1:
                sys.exit('Fatal: iw<1')
            if ie < 1:
                sys.exit('Fatal: ie<1')

            for y, mm in self.dates:
                m = (y - self.iyr1) * 12 + mm
                month = (m - 1) % 12
                kt = 0
                ts = 0.0
                for j in range(js -1, jn):
                    for i in range(iw - 1, ie):
                        if sst[i][j][m-1] < -1.77 or clim[i][j][month] == giss_data.XMISSING:
                            continue
                        kt += 1
                        ts += sst[i][j][m-1] - clim[i][j][month]
                
                box.set_value(m + moff - 1, giss_data.XMISSING)
                if kt > 0:
                    box.set_value(m + moff - 1, ts / kt)

            box.trim()
            yield box

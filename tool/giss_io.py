#!/usr/bin/env python
"""Readers and writers for GISS data.

This is just a stepping stone.
"""
__docformat__ = "restructuredtext"


import itertools
import os
import re
import struct

import math

from code import fort
import code.giss_data
import code.read_config


#: Integer code used to indicate missing data.
#:
#: This is units of 0.1 celsius. This code is only used when
#: reading or writing input and working files.
MISSING = 9999


# For all plausible integers i, float_to_tenths(tenths_to_float(i)) == i.
def tenths_to_float(t):
    if t == MISSING:
        return code.giss_data.XMISSING
    return t * 0.1


def float_to_tenths(f):
    if abs(f - code.giss_data.XMISSING) < 0.01:
        return MISSING
    return int(round(f * 10.0))


# TODO: How does this differ from float_to_tenths.
#       Answer:
#           float_to_tenths(-0.95) == -10
#           as_tenths(-0.95) == -9
def as_tenths(v):
    return int(math.floor(v * 10 + 0.5))


# TODO: Probably should be a generator.
def convert_to_tenths(celsius_series):
    """Convert a series of celsius value to integer tenths of a degree.

    :Param celsius_series:
        A list or iterable of floating point value, where each value
        represents a temperature in celsius.

    :Return:
        A new list of integer values.

    """
    t = []
    for v in celsius_series:
        if v >= code.giss_data.XMISSING:
            t.append(MISSING)
        else:
            t.append(float_to_tenths(v))
    return t


# TODO: Probably should be a generator.
def convert_tenths_to_float(tenths_series):
    """The inverse of `convert_to_tenths`."""
    t = []
    for v in tenths_series:
        if v == MISSING:
            t.append(code.giss_data.XMISSING)
        else:
            t.append(tenths_to_float(v))
    return t


def open_or_uncompress(filename):
    """Opens the text file `filename` for reading.  If this fails then
    it attempts to find a compressed version of the file by appending
    '.gz' to the name and opening that (uncompressing it on
    the fly).

    """
    try:
        return open(filename)
    except IOError:
        # When none of filename, nor filename.gz exists we
        # want to pretend that the exception comes from the original
        # call to open, above.  Otherwise the user can be confused by
        # claims that "foo.gz" does not exist when they tried to open
        # "foo".  In order to maintain this pretence, we have to get
        # the exception info and save it. See
        # http://blog.ianbicking.org/2007/09/12/re-raising-exceptions/
        import sys
        exception = sys.exc_info()
        try:
            import gzip
            return gzip.open(filename + '.gz')
        except IOError:
            pass
        raise exception[0], exception[1], exception[2]


class SubboxWriter(object):
    """Produces a GISTEMP SBBX (subbox) file; typically the output of
    step3 (and 4), and the input to step 5.
    """

    def __init__(self, rawfile, bos='>', trimmed=True):
        self.trimmed = trimmed
        self.bos = bos
        self.f = fort.File(rawfile, bos=self.bos)
        self.meta = None
        self.buf_record = None

    def add_meta(self, meta):
        self.meta = meta
        assert self.meta.mavg == 6, "Only monthly averages supported"

    def _flush(self, mo1):
        if not self.buf_record:
            m = self.meta
            rec = struct.pack(self.bos + '8i80s', mo1, m.kq, m.mavg,
                    m.monm, m.monm4, m.yrbeg, m.missing_flag,
                    m.precipitation_flag, m.title)
        else:
            b = self.buf_record
            if self.trimmed and b.station_months == 0:
                # Write as trimmed record.
                fmt = "iiiiiiif1f"
                rec = struct.pack(self.bos + fmt, mo1,
                    b.lat_S, b.lat_N, b.lon_W, b.lon_E, b.stations,
                        b.station_months, b.d, 9999.0)
            else:
                fmt = "iiiiiiif%df" % b.n
                rec = struct.pack(self.bos + fmt, mo1,
                    b.lat_S, b.lat_N, b.lon_W, b.lon_E, b.stations,
                        b.station_months, b.d, *b.series)
        self.f.writeline(rec)

    def write(self, record):
        if self.meta is None:
            assert hasattr(record, "precipitation_flag"), (
                    "First record must be SubboxMetaData")
            self.meta = record
        else:
            if self.trimmed and record.station_months == 0:
                self._flush(1)
            else:
                self._flush(record.n)
            self.buf_record = record

    def close(self):
        self._flush(1)
        self.f.close()


class StationRecordWriter(object):
    """Writes station records into a binary format that is the usual
    output of Step 2.  This format is partly documented in the header
    comment of the program STEP3/to.SBBXgrid.f from the GISTEMP sources,
    where a comment says "READS NCAR FILES".  NCAR is the National
    Center for Atmospheric Research (as US government agency).

    Step 3 takes files in this format as input.
    """
    def __init__(self, rawfile, bos='>'):
        self.bos = bos
        self.f = fort.File(rawfile, bos=self.bos)
        self.meta = None
        self.buf_record = None

    def add_meta(self, meta):
        self.meta = meta

    def _flush(self, rel_first_month, rel_last_month):
        if not self.buf_record:
            m = self.meta
            data = struct.pack(self.bos + '9i80s', rel_first_month, m.kq, m.mavg,
                    m.monm, m.monm4, m.yrbeg, m.missing_flag,
                    m.precipitation_flag, rel_last_month, m.title)

        else:
            r = self.buf_record
            s = r.station
            compound_name = "%s%c%c%c%3s" % (
                      s.name, s.US_brightness, s.pop, s.GHCN_brightness,
                      s.uid[:3])

            fmt = "%di" % len(r.series)
            data = struct.pack(self.bos + fmt, *convert_to_tenths(r.series))
            data += struct.pack(self.bos + "iiii36sii", as_tenths(s.lat),
                    as_tenths(s.lon), r.short_id, s.elevation, compound_name,
                    rel_first_month, rel_last_month)

        self.f.writeline(data)

    def write(self, record):
        if self.meta is None:
            assert hasattr(record, "precipitation_flag"), (
                    "First record must be StationMetaData")
            self.meta = record
        else:
            self._flush(record.rel_first_month, record.rel_last_month)
            self.buf_record = record

    def close(self):
        self._flush(MISSING, MISSING)
        self.f.close()


class StationReader(object):
    """Reads files in the format produced by StationRecordWriter.  Files
    is this format are the input to Step 3 (and the output from Step 2).
    """
    def __init__(self, rawfile, bos='>'):
        self.bos = bos
        self.f = fort.File(rawfile, bos=self.bos)
        rec = self.f.readline()
        (self.min_month, kq, mavg, monm, monm4, yrbeg, missing_flag,
                precipitation_flag, self.max_month,
                title) = struct.unpack(self.bos + '9i80s', rec)
        self.meta = code.giss_data.StationMetaData(self.min_month, kq, mavg,
                monm, monm4, yrbeg, missing_flag, precipitation_flag,
                self.max_month, title)

    def __iter__(self):
        return self._it()

    def _it(self):
        yield self.meta

        for rec in self.f:
            min_month, max_month = self.min_month, self.max_month
            n = max_month - min_month + 1
            fmt = "%diiiii36sii" % n
            fields = struct.unpack(self.bos + fmt, rec)
            series = fields[:n]
            (lat, lon, ident, elevation, name,
                    self.min_month, self.max_month) = fields[n:]
            country_code = name[-3:]
            uid = "%3s%09d" % (country_code, ident)
            station = code.giss_data.StationRecord(uid)
            # TODO: Rename min_month to first_month.
            # TODO: Handle magic 1880, should use meta info!
            station.set_series_from_tenths(1880 * 12 + min_month, series)

            # TODO: We could verify that the station data is correct.
            yield station


class SubboxReader(object):
    """Reads GISS subbox files (SBBX).  These files are output by Step
    3, and consumed by Step 5.  Step 4 both reads and writes a subbox
    file.
    """
    def __init__(self, rawfile, bos='>'):
        self.bos = bos
        self.f = fort.File(rawfile, bos=self.bos)
        rec = self.f.readline()
        (self.mo1, kq, mavg, monm, monm4, yrbeg, missing_flag,
                precipitation_flag,
                title) = struct.unpack(self.bos + '8i80s', rec)
        self.meta = code.giss_data.SubboxMetaData(self.mo1, kq, mavg, monm,
                monm4, yrbeg, missing_flag, precipitation_flag, title)
        assert self.meta.mavg == 6, "Only monthly averages supported"

    def info(self):
        m = self.meta
        return [self.mo1, m.kq, m.mavg, m.monm, m.monm4, m.yrbeg, m.missing_flag,
                m.precipitation_flag]

    def __iter__(self):
        return self._it()

    def _it(self):
        yield self.meta

        for rec in self.f:
            mo1 = self.mo1
            fmt = "iiiiiiif%df" % mo1
            fields = struct.unpack(self.bos + fmt, rec)
            series = fields[8:]
            (self.mo1, lat_S, lat_N, lon_W, lon_E, stations,
                    station_months, d) = fields[:8]
            subbox = code.giss_data.SubboxRecord(
                    lat_S, lat_N, lon_W, lon_E,
                    stations, station_months, d, series)
            yield subbox

    def __getattr__(self, name):
        return getattr(self.meta, name)


def StationTsReader(path):
    """Opens a file in Ts.txt file format (output by step 1) and yields
    each station.

    """

    f = open(path)

    for (record_line, lines) in itertools.groupby(
      f, lambda line: line[0] == ' '):
        if record_line: 
            # Line beginning with a blank introduces a new station
            # record.
            lines = list(lines)
            assert len(lines) == 1
            line = lines[0]
            id12 = line[10:22]
            record = code.giss_data.StationRecord(id12)
        else:
            # Lines consists of the temperature series. Year +
            # temperature values, as a set of contiguous years.
            for line in lines:
                data = [int(v) for v in line.split()]
                record.add_year(data[0], convert_tenths_to_float(data[1:]))
            yield record


class StationTsWriter(object):
    """Writes a file in Ts.txt format.  The traditional output of
    step1.
    """

    def __init__(self, path):
        self.f = open(path, "w")

    def write(self, record):
        station = record.station
        self.f.write(' %4d%5d%s%4s%-36s\n' % (
                as_tenths(station.lat), as_tenths(station.lon),
                record.uid, station.elevation, station.name))

        data = convert_to_tenths(record.series)
        for y, offset in enumerate(range(0, record.n, 12)):
            months = ["%5d" % v for v in data[offset: offset + 12]]
            self.f.write('%4d%s\n' % (y + record.first_year, ''.join(months)))

    def close(self):
        self.f.close()


def V2MeanReader(path, year_min=-9999):
    """Reads a file in GHCN v2.mean format and yields each station.
    
    If `year_min` is specified, then only years >= year_min are kept
    (the default, -9999, effectively keeps all years."""

    f = open(path)
    def id12(l):
        """Extract the 12-digit station record identifier."""
        return l[:12]

    def v2_int(s):
        v = int(s)
        if v == -9999:
            return MISSING
        return v

    for (id, lines) in itertools.groupby(f, id12):
        # lines is a set of lines which all begin with the same 12
        # character id
        record = code.giss_data.StationRecord(id)
        prev_line = None
        for line in lines:
            if line != prev_line:
                year = int(line[12:16])
                tenths = [v2_int(line[a:a+5]) for a in range(16, 16+12*5, 5)]
                if year >= year_min:
                    record.add_year(year, convert_tenths_to_float(tenths))
                prev_line = line
            else:
                print ("NOTE: repeated record found: Station %s year %s;"
                       " data are identical" % (line[:12],line[12:16]))

        if not record.is_empty():
            yield record

    f.close()


class V2MeanWriter(object):
    """Write a file in GHCN v2.mean format."""

    def __init__(self, path):
        self.f = open(path, "w")

    def to_text(self, t):
        if t == MISSING:
            return "-9999"
        else:
            return "%5d" % t

    def write(self, record):
        for year in range(record.first_year, record.last_year + 1):
            if not record.has_data_for_year(year):
                continue
            temps = [self.to_text(t)
                    for t in convert_to_tenths(record.get_a_year(year))]
            self.f.write('%s%04d%s\n' % (record.uid, year, ''.join(temps)))

    def close(self):
        self.f.close()


antarc_discard_re = re.compile(r'^$|^Get |^[12A-Z]$')
antarc_temperature_re = re.compile(r'^(.*) .* *temperature')

def read_antarctic(path, station_path, discriminator):
    stations = read_antarc_station_ids(station_path, discriminator)
    record = None
    for line in open(path):
        if antarc_discard_re.search(line):
            continue
        station_line = antarc_temperature_re.match(line)
        if station_line:
            station_name = station_line.group(1)
            station_name = station_name.replace('\\','')
            id12 = stations[station_name]
            if record is not None:
                yield record
            record = code.giss_data.StationRecord(id12)
            continue
        line = line.strip()
        if line.find('.') >= 0 and line[0] in '12':
            year, data = read_antarc_line(line)
            if year >= code.giss_data.BASE_YEAR:
                record.add_year(year, data)

    if record is not None:
        yield record


austral_discard_re = re.compile(r'^$|:')
austral_header_re = re.compile(r'^\s*(.+?)  .*(E$|E )')

def read_australia(path, station_path, discriminator):
    stations = read_antarc_station_ids(station_path, discriminator)
    record = None
    for line in open(path):
        if austral_discard_re.search(line):
            continue
        station_line = austral_header_re.match(line)
        if station_line:
            station_name = station_line.group(1).strip()
            id12 = stations[station_name]
            if record is not None:
                yield record
            record = code.giss_data.StationRecord(id12)
            continue
        line = line.strip()
        if line.find('.') >= 0 and line[0] in '12':
            year, data = read_antarc_line(line)
            if year >= code.giss_data.BASE_YEAR:
                record.add_year(year, data)

    if record is not None:
        yield record

def read_antarc_line(line):
    """Convert a single line from the Antarctic/Australasian dataset
    files into a year and a 12-tuple of floats (the temperatures in
    Centigrade).
    """

    year = int(line[:4])
    line = line[4:]
    tuple = []
    if line[6] == '.' or line[7] == '-':
        # Some of the datasets are 12f8.1 with missing values as '       -'.
        for i in range(0,12):
            tuple.append(read_float(line[i*8:i*8+8]))
    else:
        # Others are xx12f7.1 or xxx12f7.1 with missing values as '       '.
        np = line.find('.')
        if np < 0:
            raise ValueError, "Non-data line encountered: '%s'" % line
        position = (np % 7) + 2
        for i in range(0,12):
            tuple.append(read_float(line[i*7+position:i*7+7+position]))
    return (year, tuple)


def read_antarc_station_ids(path, discriminator):
    """Reads a SCARs station ID files and returns a dictionary
    mapping station name to the 12-digit station ID.
    """

    dict = {}
    for line in open(path):
        id11 = line[:11]
        station = line[12:42].strip()
        dict[station] = id11 + discriminator
    return dict


def read_USHCN(path, stations):
    # TODO: This is reading data similar to other line based records,
    #       but the algorithm is somewhat different. See examples that
    #       use itertools.groupby.
    # NOTE: Sometimes, there are two sets of data for a station.
    #       Where the same year appears in both data sets, the second year's
    #       data is used. I do not know whether this is intentional or not.
    def fill_record():
        for year, temps in sorted(years_data.iteritems()):
            record.add_year(year, temps)
        return record

    record = None
    for line in open_or_uncompress(path):
        if line[6] != '3': # 3 indicates mean temperatures
            continue
        USHCN_station = int(line[0:6])
        id12 = stations[USHCN_station]
        year = int(line[7:11])
        if year < code.giss_data.BASE_YEAR: # discard data before 1880
            continue
        if record is None or id12[:11] != record.station_uid:
            if record is not None:
                yield fill_record()
            record = code.giss_data.StationRecord(id12)
            years_data = {}

        temps = []
        valid = False
        for m in range(0,12):
            temp_fahrenheit = int(line[m*7+11:m*7+17])
            flag = line[m*7+17]
            if ((flag in 'EQ') or              # interpolated data
                (temp_fahrenheit == -9999)) :  # absent data
                temp = code.giss_data.XMISSING
            else:
                # tenths of degrees F to degrees C
                temp = (temp_fahrenheit - 320) * 5/90.0
                valid = True
            temps.append(temp)
        if valid: # some valid data found
            # We cannot add using 'add_year' because that breaks
            # things when years are duplicated in the data; (see not
            # above). The record is actually filled using ``fil_record``.
            years_data.setdefault(year, [])
            years_data[year] = temps

    if record is not None:
        yield fill_record()


def read_USHCN_stations(ushcn_v1_station_path, ushcn_v2_station_path):
    """Reads the USHCN station list and returns a dictionary
    mapping USHCN station ID to 12-digit station ID.
    """

    stations = {}
    for line in open(ushcn_v1_station_path):
        (USHCN_id, WMO_id, duplicate) = line.split()
        USHCN_id = int(USHCN_id)
        if WMO_id[0:3] != '425': # non-US country_code
            raise ValueError, "non-425 station found in ushcn.tbl: '%s'" % line
        if duplicate != '0':
            raise ValueError, "station in ushcn.tbl with non-zero duplicate: '%s'" % line
        stations[USHCN_id] = WMO_id + '0'
    # some USHCNv2 station IDs convert to USHCNv1 station IDs:
    for line in open(ushcn_v2_station_path):
        (v2_station,_,v1_station,_) = line.split()
        stations[int(v2_station)] = stations[int(v1_station)]
    return stations


def read_hohenpeissenberg(path):
    """reads the Hohenpeissenberg data from
    input/t_hohenpeissenberg_200306.txt_as_received_July17_2003
    which has a header line and then one line per year from 1781.
    We only want data from 1880 to 2002.
    """

    record = code.giss_data.StationRecord('617109620002')
    for line in open(path):
        if line[0] in '12':
            year = int(line[:4])
            if year < 1880 or year > 2002:
                continue
            data = line.split()
            temps = map(read_float, data[1:13])
            assert len(temps) == 12
            record.add_year(year, temps)
    return record


def read_float(s):
    """Returns the float converted from the argument string.
    If float conversion fails, returns XMISSING.
    """

    try:
        return float(s)
    except:
        return code.giss_data.XMISSING


# Each of the stepN_input functions below produces an iterator that
# yields that data for that step feeding from data files.
# Each of the stepN_output functions below is effectively a "tee" that
# writes the data to a file; they each take a data object (an
# iterator), write each item to a file, and yield each item.
def step0_input():
    class Struct:
        pass
    input = Struct()
    input.ushcn_stations = read_USHCN_stations('input/ushcn2.tbl', 'input/ushcnV2_cmb.tbl')
    input.ushcn_source = read_USHCN("input/9641C_201002_F52.avg", input.ushcn_stations)
    input.ghcn_source = V2MeanReader("input/v2.mean", code.giss_data.BASE_YEAR)
    input.antarc_source = itertools.chain(
            read_antarctic("input/antarc1.txt", "input/antarc1.list", '8'),
            read_antarctic("input/antarc3.txt", "input/antarc3.list", '9'),
            read_australia("input/antarc2.txt", "input/antarc2.list", '7'))
    input.hohenpeissenberg = read_hohenpeissenberg(
            "input/t_hohenpeissenberg_200306.txt_as_received_July17_2003")

    return input

def step0_output(data):
    out = V2MeanWriter("work/v2.mean_comb")
    try:
        for thing in data:
            out.write(thing)
            yield thing
    finally:
        out.close()

def step1_input():
    return V2MeanReader("work/v2.mean_comb")

def step1_output(data):
    out = StationTsWriter("work/Ts.txt")
    try:
        for thing in data:
            out.write(thing)
            yield thing
    finally:
        out.close()

def step2_input():
    return StationTsReader("work/Ts.txt")

def step2_output(data):
    out = StationRecordWriter(open("work/Ts.GHCN.CL.PA", "wb"), bos='<')
    try:
        for thing in data:
            out.write(thing)
            yield thing
    finally:
        out.close()

def step3_input():
    return StationReader(open("work/Ts.GHCN.CL.PA", "rb"), bos='<')

def step3_output(data):
    out = SubboxWriter(open('work/SBBX1880.Ts.GHCN.CL.PA.1200', 'wb'),
      trimmed=False)
    try:
        for thing in data:
            out.write(thing)
            yield thing
    finally:
        out.close()

def step4_input():
    while True:
        yield 'Step 4 dummy value'

def step4_output(data):
    # The Step 4 output is slightly unusual, it is an iterable of pairs.
    # We only want to write the records from the right-hand item (the
    # ocean data).  The left-hand items are land data, already written
    # by Step 3.
    out = SubboxWriter(open('work/SBBX.HadR2', 'wb'))
    try:
        for land,ocean in data:
            out.write(ocean)
            yield land,ocean
    finally:
        out.close()

def step5_input():
    land = SubboxReader(
      open(os.path.join('work', 'SBBX1880.Ts.GHCN.CL.PA.1200'), 'rb'))
    ocean = SubboxReader(open("work/SBBX.HadR2", "rb"))

    return itertools.izip(land, ocean)

def step5_output(data):
    # Already implicit in Step 5.
    return data

#!/usr/bin/env python
# $URL$
# $Rev$
#
# giss_io.py
#
# Paul Ollis and David Jones, 2010-03-10

"""Readers and writers for datafiles used by NASA GISS GISTEMP.

Some of these file formats are peculiar to GISS, others are defined and
used by other bodies (such as NOAA's v2.mean format).

This is just a stepping stone.
"""
__docformat__ = "restructuredtext"


import itertools
import os
import re
import struct
import glob

import math

import fort
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
        return code.giss_data.MISSING
    return t * 0.1


def float_to_tenths(f):
    if abs(f - code.giss_data.MISSING) < 0.01:
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
        if v >= code.giss_data.MISSING:
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
            t.append(code.giss_data.MISSING)
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
        if self.meta is None:
            # Not even a meta record written.  We could complain here,
            # but the cause is likely to be closing files early because
            # of some underlying problem.  Complaining here would mask
            # that problem.
            return
        if not self.buf_record:
            m = self.meta
            rec = struct.pack(self.bos + '8i80s', mo1, m.kq, m.mavg,
                    m.monm, m.monm4, m.yrbeg, m.missing_flag,
                    m.precipitation_flag, m.title)
        else:
            b = self.buf_record
            # Bounding box; to be converted to integer hundredths.
            # Conventionally the 4 elements of the box are southern
            # latitude, northern latitude, western longitude, eastern
            # longitude (but the code doesn't care).
            box = b.box
            box = [int(round(x * 100)) for x in box]
            if self.trimmed and b.station_months == 0:
                # Write as trimmed record.
                fmt = "iiiiiiif1f"
                rec = struct.pack(self.bos + fmt, mo1,
                                  box[0],
                                  box[1],
                                  box[2],
                                  box[3],
                                  b.stations, b.station_months, b.d, 9999.0)
            else:
                fmt = "iiiiiiif%df" % b.n
                rec = struct.pack(self.bos + fmt, mo1,
                                  box[0],
                                  box[1],
                                  box[2],
                                  box[3],
                                  b.stations, b.station_months, b.d, *b.series)
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


# Obsolescent.
class StationRecordWriter(object):
    """Writes station records into a binary format that is the
    traditional output of Step 2 (but more recent versions of
    ccc-gistemp use a v2.mean style output).
    
    This format is partly documented in the header comment of the
    program STEP3/to.SBBXgrid.f from the GISTEMP sources, where a
    comment says "READS NCAR FILES".  NCAR is the National Center
    for Atmospheric Research (as US government agency).

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
            if parameters.use_global_brightness:
                if s.global_brightness > 35:
                    b = '3'
                elif s.global_brightness > 10:
                    b = '2'
                else:
                    b = '1'
            else:
                b = s.US_brightness
            compound_name = "%s%c%c%c%3s" % (
                      s.name, b, s.pop, s.GHCN_brightness,
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


# Obsolescent.
class StationReader(object):
    """Reads files in the format produced by StationRecordWriter.  Files
    in this format are the traditional input to Step 3, but more modern
    ccc-gistemp uses a v2.mean style input.
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
            tenths = fields[:n]
            series = convert_tenths_to_float(tenths)
            (lat, lon, ident, elevation, name,
                    self.min_month, self.max_month) = fields[n:]
            country_code = name[-3:]
            uid = "%3s%09d" % (country_code, ident)
            station = code.giss_data.StationRecord(uid)
            # TODO: Rename min_month to first_month.
            # TODO: Handle magic 1880, should use meta info!

            station.set_series(1880 * 12 + min_month, series)

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
        """Return a length 8 sequence corresponding to the INFO array
        record in the binary file.
        """
        m = self.meta
        return [self.mo1, m.kq, m.mavg, m.monm,
                m.monm4, m.yrbeg, m.missing_flag, m.precipitation_flag]

    def __iter__(self):
        yield self.meta

        for rec in self.f:
            mo1 = self.mo1
            fmt = "iiiiiiif%df" % mo1
            fields = list(struct.unpack(self.bos + fmt, rec))
            series = fields[8:]
            self.mo1 = fields[0]
            # Make an attributes dictionary.
            # The box boundaries are fields[1:5], but we need to scale
            # them to fractional degrees first:
            for i in range(1,5):
                fields[i] /= 100.0
            attr = dict(zip(
              ['lat_S', 'lat_N', 'lon_W', 'lon_E',
              'stations', 'station_months', 'd'],
              fields[1:8]))
            attr['box'] = fields[1:5]
            subbox = code.giss_data.SubboxRecord(series,
               **attr)
            yield subbox

    def __getattr__(self, name):
        return getattr(self.meta, name)


# Obsolescent.
def StationTsReader(path):
    """Opens a file in Ts.txt file format and yields each station.
    This is the traditional output of Step 1, but modern ccc-gistemp
    uses a v2.mean style output.

    """

    f = open(path)

    def isheader(line):
        """Helper function that flags header lines."""

        return line[0] == ' '

    for (header_line, lines) in itertools.groupby(f, isheader):
        if header_line: 
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


# Obsolescent.
class StationTsWriter(object):
    """Writes a file in Ts.txt format.  The traditional output of
    step1, but modern ccc-gistemp output uses a v2.mean style output.
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
    (the default, -9999, effectively keeps all years).
    
    Traditionally this file was the output of Step 0 (and of course the
    GHCN source), but modern ccc-gistemp produces this format for the
    outputs of Steps 0, 1, and 2."""

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
    """Write a file in GHCN v2.mean format. See also V2MeanReader."""

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

def DecimalReader(path, year_min=-9999):
    """Reads a file in Decimal format and yields each station.
    
    If `year_min` is specified, then only years >= year_min are kept
    (the default, -9999, effectively keeps all years).
    """

    f = open(path)
    def id12(l):
        """Extract the 12-digit station record identifier."""
        return l[:12]

    def readt(s):
        v = float(s)
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
                if year >= year_min:
                    temps = [readt(n) for n in line[16:].split()]
                    record.add_year(year, temps)
                prev_line = line
            else:
                print ("NOTE: repeated record found: Station %s year %s;"
                       " data are identical" % (line[:12],line[12:16]))

        if not record.is_empty():
            yield record

    f.close()


class DecimalWriter(object):
    """Decimal is a novel text file format, very similar to GHCN
    v2.mean format.  Each lines conists of a 12 digit record identifier
    immediately followed by a 4 digit year, then followed by the decimal
    temperature values (in degrees C) for the 12 months of the year,
    with each temperature value preceded by a space.  Missing data are
    marked with -9999."""

    def __init__(self, path):
        self.f = open(path, "w")

    def to_text(self, t):
        if t == MISSING:
            return '-9999'
        return repr(t)

    def write(self, record):
        for year in range(record.first_year, record.last_year + 1):
            if not record.has_data_for_year(year):
                continue
            temps = [self.to_text(t) for t in record.get_a_year(year)]
            self.f.write('%s%04d %s\n' % (record.uid, year, ' '.join(temps)))

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
                temp = code.giss_data.MISSING
            else:
                # tenths of degrees F to degrees C
                temp = (temp_fahrenheit - 320) * 5/90.0
                valid = True
            temps.append(temp)
        if valid: # some valid data found
            # We cannot add using 'add_year' because that breaks
            # things when years are duplicated in the data (see note
            # above). The record is actually filled using ``fill_record``.
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
    """Reads the Hohenpeissenberg data from
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
    If float conversion fails, returns MISSING.
    """

    try:
        return float(s)
    except:
        return code.giss_data.MISSING

# Find the USHCN input file
def ushcn_input_file():
    files =  ["input/ushcnv2",
              "input/9641C_201003_F52.avg",
              "input/9641C_201002_F52.avg",
              "input/9641C_200907_F52.avg"]
    for f in files:
        if os.path.exists(f) or os.path.exists(f+'.gz'):
            return f
    raise ValueError, "no USHCN data file in input directory; run tool/preflight.py."

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
    input.ushcn_source = read_USHCN(ushcn_input_file(), input.ushcn_stations)
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
    for thing in data:
        out.write(thing)
        yield thing
    print "Step0: closing output file"
    out.close()

def step1_input():
    return V2MeanReader("work/v2.mean_comb")

def step1_output(data):
    out = V2MeanWriter("work/v2.step1.out")
    for thing in data:
        out.write(thing)
        yield thing
    print "Step1: closing output file"
    out.close()

def step2_input():
    return V2MeanReader("work/v2.step1.out")

def step2_output(data):
    out = V2MeanWriter("work/v2.step2.out")
    for thing in data:
        out.write(thing)
        yield thing
    print "Step2: closing output file"
    out.close()

def step3_input():
    return V2MeanReader("work/v2.step2.out")

def step3_output(data):
    out = SubboxWriter(open('result/SBBX1880.Ts.GHCN.CL.PA.1200', 'wb'),
      trimmed=False)
    v2out = V2MeanWriter('work/v2.step3.out')
    gotmeta = False
    for thing in data:
        out.write(thing)
        if gotmeta:
            v2out.write(thing)
        gotmeta = True
        yield thing
    print "Step3: closing output file"
    out.close()

def make_3d_array(a, b, c):
    """Create an array with three dimensions.

    Actually a list-of-lists-of-lists, but the result can be treated as an
    array with dimensions ``[a][b][c]``.

    """
    x = [0.0] * c
    arr = []
    for i in range(a):
        arr.append([list(x) for j in range(b)])

    return arr

def step4_find_monthlies(latest_year, latest_month):
    dates = {}
    filename_re = re.compile('^oiv2mon\.([0-9][0-9][0-9][0-9])([0-9][0-9])(\.gz)?$')
    for f in os.listdir('input'):
        m = filename_re.match(f)
        if m:
            year = int(m.group(1))
            month = int(m.group(2))
            if (year, month) > (latest_year, latest_month):
                if m.group(3):
                    f = f[:-3]
                dates[(year, month)] = 'input/'+f
    l = dates.items()
    l.sort()
    return l

def step4_load_sst_monthlies(latest_year, latest_month):
    files = step4_find_monthlies(latest_year, latest_month)
    if not files:
        print "No more recent sea-surface data files.\n"
        return None

    first_year = files[0][0][0]
    last_year = files[-1][0][0]
    n_years = last_year - first_year + 1

    # Read in the SST data for recent years
    sst = make_3d_array(360, 180, 12 * n_years)

    dates = []
    for (date, file) in files:
        dates.append(date)
        (year, month) = date
        f = open_or_uncompress(file)
        print "reading", file
        f = fort.File(f, bos = ">")
        f.readline() # discard first record
        data = f.readline()
        f.close()
        month = 12 * (year - first_year) + month - 1
        p = 0
        for lat in range(180):
            for long in range(360):
                v, = struct.unpack(">f", data[p:p+4])
                p += 4
                sst[long][lat][month] = v

    return sst, dates

def step4_load_clim():
    f = open_or_uncompress("input/oisstv2_mod4.clim")
    f = fort.File(f, bos='>')
    data = f.readline()
    f.close()

    clim_title = data[:80]
    clim = make_3d_array(360, 180, 12)
    p = 0
    for month in range(12):
        for lat in range(180):
            for long in range(360):
                v, = struct.unpack(">f", data[p+80:p+84])
                p += 4
                clim[long][lat][month] = v
    return clim

def step4_dummy_input():
    while True:
        yield 'Step 4 dummy value'

# This is used to extract the end month/year from the title of the SBBX file.
rTitle = re.compile(r"Monthly Sea Surface Temperature anom \(C\)"
        " Had: 1880-11/1981, oi2: 12/1981- *(\d+)/(\d+)")

def step4_input(land):
    if land is None:
        land = step4_dummy_input()
    ocean = SubboxReader(open('input/SBBX.HadR2', 'rb'))
    m = rTitle.match(ocean.meta.title)
    if m is None:
        sys.stderr.write("The title in SBBX.HadR2 does not look right\n")
        sys.stderr.write("Unable to determine end month/year from:\n")
        sys.stderr.write("  %r\n" % ocean.meta.title)
        sys.exit(1)
    end_month = int(m.group(1))
    end_year = int(m.group(2))
    monthlies = step4_load_sst_monthlies(end_year, end_month)
    return (land, ocean, monthlies)

def step4_output(data):
    # The Step 4 output is slightly unusual, it is an iterable of pairs.
    # We only want to write the records from the right-hand item (the
    # ocean data).  The left-hand items are land data, already written
    # by Step 3.
    out = SubboxWriter(open('result/SBBX.HadR2', 'wb'))
    for land,ocean in data:
        out.write(ocean)
        yield land,ocean
    print "Step4: closing output file"
    out.close()

def step5_input():
    land = SubboxReader(open('result/SBBX1880.Ts.GHCN.CL.PA.1200', 'rb'))
    ocean = SubboxReader(open('result/SBBX.HadR2', 'rb'))

    return itertools.izip(land, ocean)

def step5_bx_output(data):
    bos = '>'
    box = open('result/BX.Ts.ho2.GHCN.CL.PA.1200', 'wb')
    boxf = fort.File(box, bos=bos)
    (info, title) = data.next()
    boxf.writeline(struct.pack('%s8i' % bos, *info) + title)
    yield (info, title)

    for record in data:
        (avgr, wtr, ngood, box) = record
        n = len(avgr)
        fmt = '%s%df' % (bos, n)
        boxf.writeline(struct.pack(fmt, *avgr) +
                       struct.pack(fmt, *wtr) +
                       struct.pack('%si' % bos, ngood) +
                       struct.pack('%s4i' % bos, *box))
        yield record
    print "Step5: Closing box file"
    boxf.close()

def step5_zone_titles():
    # Boundaries (degrees latitude, +ve North) of the 8 basic belts.
    band = ['90.0 N',
            '64.2 N',
            '44.4 N',
            '23.6 N',
            'EQUATOR',
            '23.6 S',
            '44.4 S',
            '64.2 S',
            '90.0 S']
    # Accumulate the titles here.
    titles = ['  LATITUDE BELT FROM %7s TO %7s' % (band[j+1], band[j]) for j in range(8)]

    titles += [
        '  NORTHERN LATITUDES: 23.6 N TO  90.0 N',
        '       LOW LATITUDES: 23.6 S TO  23.6 N',
        '  SOUTHERN LATITUDES: 90.0 S TO  23.6 S',
        'NORTHERN HEMISPHERE',
        'SOUTHERN HEMISPHERE',
        'GLOBAL MEANS']

    # Ensure all titles are 80 characters long.
    titles = map(lambda s: ('%-80s' % s)[:80], titles)
    return titles


def step5_output(data):
    (info, data, wt, ann, monmin, title) = data
    XBAD = 9999
    iy1tab = 1880

    zone_titles = step5_zone_titles()

    titl2 = ' zones:  90->64.2->44.4->23.6->0->-23.6->-44.4->-64.2->-90                      '
    iyrbeg = info[5]
    jzm = len(ann)
    iyrs = len(ann[0])
    monm = iyrs * 12
    
    out = ['ZonAnn', 'GLB', 'NH', 'SH']
    out = [open('result/'+bit+'.Ts.ho2.GHCN.CL.PA.txt', 'w')
            for bit in out]
    zono = open('result/ZON.Ts.ho2.GHCN.CL.PA.1200', 'wb')

    bos = '>'

    zono = fort.File(zono, bos)

    # Create and write out the header record of the output files.
    print >> out[0], ' Annual Temperature Anomalies (.01 C) - ' + title[28:80]
    for f in out[1:]:
        print >> f, title

    # iord literal borrowed exactly from Fortran...
    iord = [14,12,13, 9,10,11, 1,2,3,4,5,6,7,8]
    # ... and then adjusted for Python index convention.
    iord = map(lambda x: x-1, iord)
    # Display the annual means.
    def annasstr(z):
        """Helper function that returns the annual anomaly for zone *z*
        as a string representation of an integer (the integer is the
        anomaly scaled by 100 to convert to centikelvin).

        The returned value is a string that is 5 characters long.  If
        the integer will not fit into a 5 character string, '*****' is
        returned (this emulates the Fortran convention of formatting
        999900 (which is the XBAD value in centikelvin) as a '*****'.
        
        The year, *iy*, is lexically captured which is a bit horrible.
        """
        x = int(math.floor(100*ann[z][iy] + 0.5))
        x = '%5d' % x
        if len(x) > 5:
            return '*****'
        return x

    iyrsp = iyrs
    # Check (and skip) incomplete year.
    if data[-1][-1][-1] > 8000:
        iyrsp -= 1
    banner = """
                           24N   24S   90S     64N   44N   24N   EQU   24S   44S   64S   90S
Year  Glob  NHem  SHem    -90N  -24N  -24S    -90N  -64N  -44N  -24N  -EQU  -24S  -44S  -64S Year
""".strip('\n')
    for iy in range(iy1tab - iyrbeg, iyrsp):
        if (iy+iyrbeg >= iy1tab+5 and ((iy+iyrbeg) % 20 == 1) or
          iy == iy1tab - iyrbeg):
            print >> out[0]
            print >> out[0], banner
        iyr = iyrbeg+iy
        print >> out[0], ('%4d' + ' %s'*3 + '  ' + ' %s'*3 +
                          '  ' + ' %s'*8 + '%5d') % tuple([iyr] +
          [annasstr(iord[zone]) for zone in range(jzm)] + [iyr])
    # The trailing banner is just like the repeated banner, except that
    # "Year  Glob  NHem  SHem" appears on on the first line, not the
    # second line (and the same for the "Year" that appears at the end
    # of the line).  *sigh*.
    banner = banner.split('\n')
    banner[0] = banner[1][:24] + banner[0][24:] + ' Year'
    banner[1] = ' '*24 + banner[1][24:-5]
    banner = '\n'.join(banner)
    print >> out[0], banner
    print >> out[0]

    tit = ['    GLOBAL','N.HEMISPH.','S.HEMISPH.']
    # Shift the remaining 3 output files so that the indexing works out.
    out = out[1:]
    banner = 'Year   Jan  Feb  Mar  Apr  May  Jun  Jul  Aug' + \
      '  Sep  Oct  Nov  Dec    J-D D-N    DJF  MAM  JJA  SON  Year'
    # All the "WRITE(96+J" stuff in the Fortran is replaced with this
    # enumeration into the *out* array (an array of file descriptors).
    for j,outf in enumerate(out):
        print >> outf, (tit[j] + ' Temperature Anomalies' + 
          ' in .01 C     base period: 1951-1980')
        for iy in range(iy1tab-iyrbeg, iyrs):
            iout = [100*XBAD]*18
            if (iy+iyrbeg >= iy1tab+5 and ((iy+iyrbeg) % 20 == 1) or
              iy == iy1tab - iyrbeg):
                print >> outf
                print >> outf, banner
            # *data* for this zone, avoids some duplication of code.
            zdata = data[iord[j]]
            # :todo: Would probably be better to have a little 4-long
            # seasonal array to do the computation in.
            awin = 9999
            if iy > 0:
                awin = zdata[iy-1][11] + zdata[iy][0] + zdata[iy][1]
            aspr = sum(zdata[iy][2:5])
            asmr = sum(zdata[iy][5:8])
            afl  = sum(zdata[iy][8:11])
            if awin < 8000:
                iout[14] = int(round(100.0*awin/3))
            if aspr < 8000:
                iout[15] = int(round(100.0*aspr/3))
            if asmr < 8000:
                iout[16] = int(round(100.0*asmr/3))
            if afl < 8000:
                iout[17] = int(round(100.0*afl/3))
            ann2=awin+aspr+asmr+afl
            if ann2 < 8000:
                iout[13] = int(round(100.0*ann2/12))
            ann1=ann[iord[j]][iy]
            if iy == iyrs-1 and zdata[iy][-1] > 8000:
                ann1 = 9999
            if ann1 < 8000:
                iout[12] = int(round(100.0*ann[iord[j]][iy]))
            for m in range(12):
                iout[m] = int(round(100.0*zdata[iy][m]))
            iyr = iyrbeg+iy
            # Convert each of *iout* to a string, storing the results in
            # *sout*.
            sout = [None]*len(iout)
            for i,x in enumerate(iout):
                # All the elements of iout are formatted to width 5,
                # except for element 13 (14 in the Fortran code), which
                # is length 4.
                if i == 13:
                    x = '%4d' % x
                    if len(x) > 4:
                        x = '****'
                else:
                    x = '%5d' % x
                    if len(x) > 5:
                        x = '*****'
                sout[i] = x
            print >> outf, (
              '%4d ' + '%s'*12 + '  %s%s  ' + '%s'*4 + '%6d') % tuple(
              [iyr] + sout + [iyr])
        print >> outf, banner

    # Save monthly means on disk.
    zono.writeline(struct.pack(bos + '8i', *info) +
                   title + titl2)

    fmt_mon = bos + '%df' % monm
    for jz in range(jzm):
        zono.writeline(struct.pack(fmt_mon, *itertools.chain(*data[jz])) +
                       struct.pack(fmt_mon, *itertools.chain(*wt[jz])) +
                       zone_titles[jz])
    return "Step 5 Completed"

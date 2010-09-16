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


import glob
import itertools
import math
import os
import re
import struct


# Clear Climate Code
import extend_path
import fort
import code.giss_data
import code.read_config


#: Integer code used to indicate missing data.
#:
#: This is units of 0.1 celsius. This code is only used when
#: reading or writing input and working files.
MISSING = 9999

# For all plausible integers, converting them from external to internal
# to external again will preserve their value.
def tenths_to_float(t):
    if t == MISSING:
        return code.giss_data.MISSING
    return t * 0.1

# TODO: How does this differ from float_to_tenths.
#       Answer:
#           float_to_tenths(-0.95) == -10
#           as_tenths(-0.95) == -9
# Note: only used by obsolete file formats.
def as_tenths(v):
    return int(math.floor(v * 10 + 0.5))

# Load GHCN V2 metadata, if present.
v2inv = os.path.join('input', 'v2.inv')
try:
    v2meta = code.giss_data.read_stations(v2inv, format='v2')
except IOError:
    v2meta = None


# TODO: Probably should be a generator.
def internal_to_external(series, scale=0.1):
    """Convert a series of values to external representation by
    converting to integer tenths (or other scale).  Normally
    this is used to convert a series from degrees Celcius to tenths
    of a degree.

    :Param series:
        A list or iterable of floating point value; usually each value
        represents a temperature in Celsius.

    :Return:
        A new list of values (ints).

    """

    # Note: 1/0.1 == 10.0; 1/0.01 == 100.0 (in other words even though
    # 0.1 and 0.01 are not stored exactly, their reciprocal is exactly
    # an integer)
    scale = 1.0/scale

    def toint(f):
        # :todo: Use of abs() probably not needed.
        if abs(f - code.giss_data.MISSING) < 0.01:
            return MISSING
        return int(round(f * scale))

    return [toint(v) for v in series]

# TODO: Probably should be a generator.
def convert_tenths_to_float(tenths_series):
    """The inverse of `internal_to_external`."""
    return [tenths_to_float(v) for v in tenths_series]


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
                fmt = "iiiiiiif%df" % len(b)
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
                self._flush(len(record))
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

            # The 11-digit station identifier with the 3-digit country
            # code at the front removed.  Only used in this output
            # format.
            short_id = r.uid[3:]

            fmt = "%di" % len(r.series)
            data = struct.pack(self.bos + fmt, *internal_to_external(r.series))
            data += struct.pack(self.bos + "iiii36sii", as_tenths(s.lat),
                    as_tenths(s.lon), int(short_id), s.elevation,
                    compound_name, rel_first_month, rel_last_month)

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
            station = code.giss_data.Series(uid=uid)
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
            subbox = code.giss_data.Series(series=series, **attr)
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
            record = code.giss_data.Series(uid=id12)
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

        data = internal_to_external(record.series)
        for y, offset in enumerate(range(0, len(record), 12)):
            months = ["%5d" % v for v in data[offset: offset + 12]]
            self.f.write('%4d%s\n' % (y + record.first_year, ''.join(months)))

    def close(self):
        self.f.close()

def GHCNV3Reader(inp, meta=None, year_min=None, scale=None):
    """Reads a file in GHCN V3 .dat format and yields each station
    record (as a giss_data.Series instance).  For now, this treats
    all the data for a station as a single record (contrast with GHCN V2
    which could have several "duplicates" for a single station).

    If a *meta* dict is supplied then the Series instance will have its
    "station" attribute set to value corresponding to the 11-digit ID in
    the *meta* dict.

    If `year_min` is specified, then only years >= year_min are kept
    (the default, None, keeps all years).
    
    If *scale* is specified then the (integer) values in the file are
    multiplied by *scale* before being returned.  When it is not
    specified (the normal case), the scale is derived from the element
    specified in the file (normally for monthly means this is "TAVG" and
    the scale implied by that is 0.01 (degrees C)).

    See ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/v3/README.pdf for format
    of this file.
    """

    def id11(l):
        """Extract the 11-digit station identifier."""
        return l[:11]

    noted_element = False
    def note_element(element):
        """Print the meteorological element we are reading."""
        friendly = dict(TAVG='average temperature')
        print "(Reading %s)" % friendly[element]

    element_scale = dict(TAVG=0.01)
    # Flags that cause value to be rejected. :todo: make parameter.
    reject = 'DKOSTW'

    def convert(s):
        """Convert single value. *s* is the 8 character string: 5
        characters for value, 3 for flags."""

        # This function captures *multiplier* which can, in principle,
        # change for each line.

        v = int(s[:5])
        # Flags for Measurement (missing days), Quality, and
        # Source.
        m,q,s = s[5:8]
        if q in reject or v == -9999:
            v = MISSING
        else:
            v *= multiplier
        return v

    all_missing = [MISSING]*12

    for id,lines in itertools.groupby(inp, id11):
        key = dict(uid=id+'G',
          first_year=year_min,
          )
        if meta and meta.get(id):
            key['station'] = meta[id]
        record = code.giss_data.Series(**key)
        for line in lines:
            year = int(line[11:15])
            element = line[15:19]
            if not noted_element:
                note_element(element)
                noted_element = True
            if scale:
                multiplier = scale
            else:
                multiplier = element_scale[element]
            values = [convert(line[i:i+8]) for i in range(19,115,8)]
            if values != all_missing:
                record.add_year(year, values)
        if len(record) != 0:
            yield record


def V2MeanReader(path=None, file=None, meta=None, year_min=None):
    """Reads a file in GHCN v2.mean format and yields each station.

    If a *meta* dict is supplied then the Series instance will have its
    "station" attribute set to value corresponding to the 11-digit ID in
    the *meta* dict.
    
    If `year_min` is specified, then only years >= year_min are kept
    (the default, None, keeps all years).

    Traditionally a file in this format was the output of Step 0 (and
    of course the format used by the GHCN source), but modern ccc-gistemp
    produces this format for the outputs of Steps 0, 1, and 2."""

    if path:
        f = open(path)
    else:
        f = file

    def id12(l):
        """Extract the 12-digit station record identifier."""
        return l[:12]

    def v2_float(s):
        """Convert a single datum from string to float; converts missing
        values from their V2 representation, "-9999", to internal
        representation, giss_data.MISSING; scales temperatures to
        convert them from integer tenths to fractional degrees C.
        """

        if "-9999" == s:
            return code.giss_data.MISSING
        else:
            return float(s) * 0.1

    # The Series.add_year protocol assumes that the years are in
    # increasing sequence.  This is so in the v2.mean file but does not
    # seem to be documented (it seems unlikely to change either).

    # Group the input file into blocks of lines, all of which share the
    # same 12-digit ID.
    for (id, lines) in itertools.groupby(f, id12):
        key = dict(uid=id, first_year=year_min)
        # 11-digit station ID.
        stid = id[:11]
        if meta and meta.get(stid):
            key['station'] = meta[stid]
        record = code.giss_data.Series(**key)
        prev_line = None
        for line in lines:
            if line != prev_line:
                year = int(line[12:16])
                temps = [v2_float(line[a:a+5]) for a in range(16, 16+12*5, 5)]
                record.add_year(year, temps)
                prev_line = line
            else:
                print ("NOTE: repeated record found: Station %s year %s;"
                       " data are identical" % (line[:12],line[12:16]))

        if len(record) != 0:
            yield record

    f.close()


class V2MeanWriter(object):
    """Write a file in GHCN v2.mean format. See also V2MeanReader."""

    def __init__(self, path=None, file=None, scale=0.1, **k):
        if path is not None:
            self.f = open(path, "w")
        else:
            self.f = file
        self.scale = scale

    def to_text(self, t):
        if t == MISSING:
            return "-9999"
        else:
            return "%5d" % t

    def write(self, record):
        """Write an entire record out."""
        for year in range(record.first_year, record.last_year + 1):
            if not record.has_data_for_year(year):
                continue
            self.writeyear(record.uid, year, record.get_a_year(year))

    def writeyear(self, uid, year, temps):
        """Write a single year's worth of data out.  *temps* should
        contain 12 monthly values."""

        strings = [self.to_text(t)
                   for t in internal_to_external(temps, scale=self.scale)]
        self.f.write('%s%04d%s\n' % (uid, year, ''.join(strings)))

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
        record = code.giss_data.Series(uid=id)
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

        if len(record) != 0:
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


# :todo: read_antarctic and read_australia are probably too similar.
# :todo: read_antarctic and read_australia would probably benefit from
# itertools.groupby

antarc_discard_re = re.compile(r'^$|^Get |^[12A-Z]$')
antarc_temperature_re = re.compile(r'^(.*) .* *temperature')

def read_antarctic(path, station_path, discriminator, meta=None):
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
            key = dict(uid=id12)
            id11 = id12[:11]
            if meta and meta.get(id11):
                key['station'] = meta[id11]
            record = code.giss_data.Series(**key)
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

def read_australia(path, station_path, discriminator, meta=None):
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
            key = dict(uid=id12)
            id11 = id12[:11]
            if meta and meta.get(id11):
                key['station'] = meta[id11]
            record = code.giss_data.Series(**key)
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


def convert_USHCN_id(record_stream, stations):
    """Convert all the records in *record_stream* from having (6-digit)
    USHCN identifiers to having 12-digit GHCN identifiers (using the
    *stations* dictionary)."""

    for record in record_stream:
        id12 = stations[int(record.uid)]
        record.uid = id12
        yield record

def convert_F_to_C(record_stream):
    """Convert each of the series from degrees Fahrenheit to degrees
    Celsius."""

    def convert_datum(x):
        if code.giss_data.invalid(x):
            return x
        # degrees F to degrees C
        return (x - 32) * 5 / 9.0

    for record in record_stream:
        record.set_series(record.first_month,
          map(convert_datum, record.series))
        yield record

def read_USHCN(path):
    """Open the USHCN V2 file *path* and yield a series of temperature
    records.  Each record is in degrees Fahrenheit (the unit used in the
    USHCN files), and will have its `uid` attribute set to its USHCN
    identifier.

    Data marked as missing (-9999 in the USHCN file) or flagged with 'E'
    or 'Q' will be replaced with MISSING.
    """

    def id6(l):
        """The 6 digit USHCN identifier."""
        return l[0:6]

    noted_element = False
    def note_element(line):
        """Print the meteorological element we are reading."""
        # See ftp://ftp.ncdc.noaa.gov/pub/data/ushcn/v2/monthly/readme.txt
        element = line[6]
        assert element in '1234'
        element = {'1':'mean maximum temperature',
                   '2':'mean minimum temperature',
                   '3':'average temperature',
                   '4':'precipitation',
                  }[element]
        print "(Reading %s)" % element

    for id,lines in itertools.groupby(open_or_uncompress(path), id6):
        record = code.giss_data.Series(uid=id)
        for line in lines:
            if not noted_element:
                note_element(line)
                noted_element = True
            # '1', '2', '3' indicate (max, min, average) temperatures.
            assert line[6] in '123'
            year = int(line[7:11])
            if year < code.giss_data.BASE_YEAR: # discard data before 1880
                continue

            temps = []
            valid = False
            for m in range(0,12):
                temp_fahrenheit = int(line[m*7+11:m*7+17])
                flag = line[m*7+17]
                if ((flag in 'EQ') or              # interpolated data
                    (temp_fahrenheit == -9999)):   # absent data
                    temp = code.giss_data.MISSING
                else:
                    # Convert to (fractional) degrees F
                    temp = temp_fahrenheit / 10.0
                    valid = True
                temps.append(temp)
            if valid: # some valid data found
                record.add_year(year, temps)
        yield record

def read_USHCN_converted(path, stations):
    """Read the USHCN data in the file *path*, converting each record to
    degrees Celsius, and converting their station identifiers to use the
    12-digit GHCN identifiers specified in the *stations* dict.
    """

    ushcn = read_USHCN(path)
    celsius = convert_F_to_C(ushcn)
    ghcn_ids = convert_USHCN_id(celsius, stations)
    return ghcn_ids


def read_USHCN_stations(ushcn_v1_station_path, ushcn_v2_station_path):
    """Reads the USHCN station list and returns a dictionary
    mapping USHCN station ID (integer) to 12-digit GHCN record ID
    (string).
    """

    stations = {}
    for line in open(ushcn_v1_station_path):
        (USHCN_id, id11, duplicate) = line.split()
        USHCN_id = int(USHCN_id)
        if not id11.startswith('425'):
            # non-US country_code
            raise ValueError, "non-425 station found in ushcn.tbl: '%s'" % line
        if duplicate != '0':
            raise ValueError, "station in ushcn.tbl with non-zero duplicate: '%s'" % line
        stations[USHCN_id] = id11 + '0'
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

    record = code.giss_data.Series(uid='617109620002')
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

def ushcn_input_file():
    """Find the USHCN input file."""
    files =  ["input/ushcnv2",
              "input/9641C_201003_F52.avg",
              "input/9641C_201002_F52.avg",
              "input/9641C_200907_F52.avg"]
    for f in files:
        if os.path.exists(f) or os.path.exists(f+'.gz'):
            return f
    raise ValueError, "no USHCN data file in input directory; run tool/preflight.py."

def ghcn3_input_file():
    """Find the GHCN V3 input file.  Looks in the input/ directory for
    sub-directories that start "ghcnm", then picks the most recent one
    and looks inside that for the .dat file (using the last 8 characters
    of the directory name to determine most recent)."""

    dirs = filter(os.path.isdir, glob.glob('input/ghcnm*'))
    dirs.sort(key=lambda x: x[:8], reverse=True)
    if not dirs:
        return None
    dir = dirs[0]
    dats = glob.glob(os.path.join(dir, '*.dat'))
    return dats[0]


# Each of the stepN_input functions below produces an iterator that
# yields that data for that step feeding from data files.
# Each of the stepN_output functions below is effectively a "tee" that
# writes the data to a file; they each take a data object (an
# iterator), write each item to a file, and yield each item.
def step0_input():
    class Struct:
        pass
    input = Struct()
    ushcn_map = read_USHCN_stations('input/ushcn2.tbl',
      'input/ushcnV2_cmb.tbl')
    input.ushcn = read_USHCN_converted(ushcn_input_file(),
      ushcn_map)
    input.ghcn = V2MeanReader("input/v2.mean",
        meta=v2meta,
        year_min=code.giss_data.BASE_YEAR)
    ghcn3file = ghcn3_input_file()
    if ghcn3file:
        invfile = ghcn3file.replace('.dat', '.inv')
        input.v3 = GHCNV3Reader(open(ghcn3file),
          meta=code.giss_data.read_stations(invfile, format='v3'),
          year_min=code.giss_data.BASE_YEAR)
    input.scar = itertools.chain(
        read_antarctic("input/antarc1.txt", "input/antarc1.list", '8',
          meta=v2meta),
        read_antarctic("input/antarc3.txt", "input/antarc3.list", '9',
          meta=v2meta),
        read_australia("input/antarc2.txt", "input/antarc2.list", '7',
          meta=v2meta))
    input.hohenpeissenberg = read_hohenpeissenberg(
            "input/t_hohenpeissenberg_200306.txt_as_received_July17_2003")

    return input

def step0_output(data):
    out = V2MeanWriter(path="work/v2.mean_comb")
    for thing in data:
        out.write(thing)
        yield thing
    print "Step0: closing output file"
    out.close()

def step1_input():
    return V2MeanReader("work/v2.mean_comb", meta=v2meta)

def step1_output(data):
    out = V2MeanWriter(path="work/v2.step1.out")
    for thing in data:
        out.write(thing)
        yield thing
    print "Step1: closing output file"
    out.close()

def step2_input():
    return V2MeanReader("work/v2.step1.out", meta=v2meta)

def step2_output(data):
    out = V2MeanWriter(path="work/v2.step2.out")
    for thing in data:
        out.write(thing)
        yield thing
    print "Step2: closing output file"
    out.close()

def step3_input():
    return V2MeanReader("work/v2.step2.out", meta=v2meta)

STEP3_OUT = os.path.join('result', 'SBBX1880.Ts.GHCN.CL.PA.1200')

def step3_output(data):
    out = SubboxWriter(open(STEP3_OUT, 'wb'),
      trimmed=False)
    v2out = V2MeanWriter(path='work/v2.step3.out', scale=0.01)
    gotmeta = False
    for thing in data:
        out.write(thing)
        if gotmeta:
            v2out.write(thing)
        gotmeta = True
        yield thing
    print "Step3: closing output file"
    out.close()
    v2out.close()

def step3c_input():
    """Use the output from the ordinary Step 3."""

    land = SubboxReader(open(STEP3_OUT, 'rb'))
    return iter(land)

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
    # Dummy metadata first.  step5.py recognises this value.
    yield None
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

def step5_input(data):
    if not data:
        land = SubboxReader(open(STEP3_OUT, 'rb'))
        ocean = SubboxReader(open('result/SBBX.HadR2', 'rb'))
        data = itertools.izip(land, ocean)
    else:
        data = ensure_landocean(data)
    # Add optional mask.
    try:
        p = os.path.join('input', 'step5mask')
        mask = open(p)
        print "Using mask from", p
    except IOError:
        mask = None
    meta = data.next()
    if mask is None:
        yield (None,) + tuple(meta)
        for land, ocean in data:
            yield None, land, ocean
    else:
        yield ('mask from %s' % mask.name,) + tuple(meta)
        for maskrow, (land, ocean) in itertools.izip(mask, data):
            maskv = float(maskrow[16:21])
            yield maskv, land, ocean

def ensure_landocean(data):
    """Ensure that the data stream has a land and an ocean series.  If
    we are piping Step 3 straight into Step 5 then we only have a land
    series.  In that case we synthesize missing ocean data.
    """

    # First item from iterator is normally a pair of metadata objects,
    # one for land, one for ocean.  If we are piping step3 straight into
    # step5 then it is not a pair.

    meta = data.next()
    try:
        land_meta, ocean_meta = meta
    except (TypeError, ValueError):
        # Use the land meta object for both land and ocean data
        land_meta,ocean_meta = meta, meta
        print "No ocean data; using land data only"
        data = add_blank(data, 'ocean')

    if land_meta is None:
        # Synthesize land data
        land_meta = ocean_meta
        print "No land data; using ocean data only"
        data = add_blank(extract_ocean(data), 'land')

    yield land_meta, ocean_meta
    for series in data:
        yield series

def extract_ocean(data):
    for land,ocean in data:
        yield ocean

def add_blank(data, required):
    """Augment a single data series with blank data to make a data
    series pair.  *required* should be 'land' to synthesize the first of
    the pair; or 'ocean' to synthesize the second of the pair.
    """

    assert required in ('land', 'ocean')

    for this_box in data:
        other_series = [MISSING] * len(this_box.series)
        other_box = code.giss_data.Series(series=other_series,
            box=this_box.box,
            stations=0, station_months=0,
            d=MISSING)
        if required == 'land':
            yield other_box, this_box
        else:
            yield this_box, other_box


def step5_bx_output(data):
    """Write box (BX) output file."""

    bos = '>'
    box = open('result/BX.Ts.ho2.GHCN.CL.PA.1200', 'wb')
    boxf = fort.File(box, bos=bos)
    info, title = data.next()
    boxf.writeline(struct.pack('%s8i' % bos, *info) + title)
    yield (info, title)

    for record in data:
        avgr, wtr, ngood, box = record
        n = len(avgr)
        fmt = '%s%df' % (bos, n)
        boxf.writeline(struct.pack(fmt, *avgr) +
                       struct.pack(fmt, *wtr) +
                       struct.pack('%si' % bos, ngood) +
                       struct.pack('%s4i' % bos, *box))
        yield record
    print "Step5: Closing box file"
    boxf.close()

def step5_mask_output(data):
    """Output the landmask used by Step 5."""

    # metadata
    yield data.next()

    out = open(os.path.join('work', 'step5mask'), 'w')

    for datum in data:
        mask,land,ocean = datum
        assert boxid_eq(land.uid, ocean.uid)
        out.write("%sMASK%.3f\n" % (land.uid, mask))
        yield datum
    out.close()

def boxid_eq(uid1, uid2):
    """Compare two IDs both of which are from subboxes.  They should be
    of the form -SS.S+EEE.EC.  They should be the same, but due to
    rounding differences, they can differ by 1 in the last digit."""

    lat1 = float(uid1[:5])
    lon1 = float(uid1[5:11])
    lat2 = float(uid2[:5])
    lon2 = float(uid2[5:11])
    return abs(lat1 - lat2) < 0.15 and abs(lon1 - lon2) < 0.15

def step5_zone_titles():
    """Return the titles used for the 14 zones."""

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
    titles = ['  LATITUDE BELT FROM %7s TO %7s' % (band[j+1], band[j])
      for j in range(8)]

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
    """Generate final Step 5 output files."""

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

    bos = '>'

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
            # Each year is formatted as a row of 18 numbers (12 months,
            # 2 different annual anomalies, and 4 seasonal).
            row = [100*XBAD]*18
            if (iy+iyrbeg >= iy1tab+5 and ((iy+iyrbeg) % 20 == 1) or
              iy == iy1tab - iyrbeg):
                print >> outf
                print >> outf, banner
            # *data* for this zone, avoids some duplication of code.
            zdata = data[iord[j]]
            # 4 seasons.
            season = [9999]*4
            if iy > 0:
                season[0] = zdata[iy-1][11] + zdata[iy][0] + zdata[iy][1]
            for s in range(1, 4):
                season[s] = sum(zdata[iy][s*3-1:s*3+2])
            # Each season slots into slots 14 to 17 of *row*.
            for s,x in enumerate(season):
                if x < 8000:
                    row[14+s] = int(round(100.0*x/3))
            # Meteorological Year is average of 4 seasons.
            metann = sum(season)
            if metann < 8000:
                row[13] = int(round(100.0*metann/12))
            # Calendar year as previously computed.
            calann = ann[iord[j]][iy]
            # For final year of data, suppress annual anomaly value unless
            # December is present (assuming a full year of data?).
            if iy == iyrs-1 and zdata[iy][-1] > 8000:
                calann = 9999
            if calann < 8000:
                row[12] = int(round(100.0*ann[iord[j]][iy]))
            # Fill in the months.
            for m in range(12):
                row[m] = int(round(100.0*zdata[iy][m]))
            # Convert each of *row* to a string, storing the results in
            # *sout*.
            formatted_row = [None]*len(row)
            for i,x in enumerate(row):
                # All the elements of *row* are formatted to width 5,
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
                formatted_row[i] = x
            year = iyrbeg+iy
            print >> outf, (
              '%4d ' + '%s'*12 + '  %s%s  ' + '%s'*4 + '%6d') % tuple(
              [year] + formatted_row + [year])
        print >> outf, banner

    # Save monthly means on disk.
    zono = open('result/ZON.Ts.ho2.GHCN.CL.PA.1200', 'wb')
    zono = fort.File(zono, bos)
    zono.writeline(struct.pack(bos + '8i', *info) +
                   title + titl2)

    fmt_mon = bos + '%df' % monm
    for jz in range(jzm):
        zono.writeline(struct.pack(fmt_mon, *itertools.chain(*data[jz])) +
                       struct.pack(fmt_mon, *itertools.chain(*wt[jz])) +
                       zone_titles[jz])
    zono.close()
    return "Step 5 Completed"

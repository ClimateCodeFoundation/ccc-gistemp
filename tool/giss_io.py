#!/usr/bin/env python
"""Readers and writers for GISS data.

This is just a stepping stone.
"""
__docformat__ = "restructuredtext"


import struct
import itertools
import re

import math

from code import fort
import code.giss_data
import code.read_config


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

    def add_record(self, record):
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

            fmt = "%di" % len(r.series_as_tenths)
            data = struct.pack(self.bos + fmt, *r.series_as_tenths)
            data += struct.pack(self.bos + "iiii36sii", s.lat_as_tenths,
                    s.lon_as_tenths, r.short_id, s.elevation, compound_name,
                    rel_first_month, rel_last_month)

        self.f.writeline(data)

    def add_record(self, record):
        if self.meta is None:
            assert hasattr(record, "precipitation_flag"), (
                    "First record must be StationMetaData")
            self.meta = record
        else:
            self._flush(record.rel_first_month, record.rel_last_month)
            self.buf_record = record

    def close(self):
        self._flush(code.giss_data.MISSING, code.giss_data.MISSING)
        self.f.close()


class StationReader(object):
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
                record.add_year_of_tenths(data[0], data[1:])
            yield record


class StationTsWriter(object):
    def __init__(self, path):
        self.f = open(path, "w")

    def add_record(self, record):
        station = record.station
        self.f.write(' %4d%5d%s%4s%-36s\n' % (
                station.lat_as_tenths, station.lon_as_tenths,
                record.uid, station.elevation, station.name))

        data = record.series_as_tenths
        for y, offset in enumerate(range(0, record.n, 12)):
            months = ["%5d" % v for v in data[offset: offset + 12]]
            self.f.write('%4d%s\n' % (y + record.first_year, ''.join(months)))

    def close(self):
        self.f.close()


def V2MeanReader(path):
    """Reads a file in GHCN v2.mean format and yields each station."""

    f = open(path)
    def id12(l):
        """Extract the 12-digit station record identifier."""
        return l[:12]

    for (id, lines) in itertools.groupby(f, id12):
        # lines is a set of lines which all begin with the same 12
        # character id
        record = code.giss_data.StationRecord(id)
        prev_line = None
        for line in lines:
            if line != prev_line:
                year = int(line[12:16])
                tenths = [int(line[a:a+5]) for a in range(16, 16+12*5, 5)]
                if year >= code.giss_data.BASE_YEAR:
                    record.add_year_of_tenths(year, tenths)
                prev_line = line
            else:
                print ("NOTE: repeated record found: Station %s year %s;"
                       " data are identical" % (line[:12],line[12:16]))

        if not record.is_empty():
            yield record

    f.close()


class V2MeanWriter(object):
    def __init__(self, path):
        self.f = open(path, "w")

    def to_text(self, t):
        if t == code.giss_data.MISSING:
            return "-9999"
        else:
            return "%5d" % t

    def add_record(self, record):
        s = self.to_text
        for year in range(record.first_year, record.last_year + 1):
            if not record.has_data_for_year(year):
                continue
            temps = [s(t) for t in record.get_a_year_as_tenths(year)]
            self.f.write('%s%04d%s\n' % (record.uid, year, ''.join(temps)))

    def close(self):
        self.f.close()


class AntarcticReader(object):
    def __init__(self, path, station_path, discriminator):
        self.f = open(path)
        self.stations = read_antarc_station_ids(station_path, discriminator)

    def __iter__(self):
        return self._it()

    def _it(self):
        record = None
        stations = self.stations
        for line in self.f:
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
                    record.add_year_of_tenths(year, data)

        if record is not None:
            yield record
        self.f.close()


class AustroAntarcticReader(AntarcticReader):
    def _it(self):
        # TODO: We can make the austro and other reader much more common.
        record = None
        stations = self.stations
        for line in self.f:
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
                    record.add_year_of_tenths(year, data)

        if record is not None:
            yield record
        self.f.close()


austral_discard_re = re.compile(r'^$|:')
austral_header_re = re.compile(r'^\s*(.+?)  .*(E$|E )')


def round_to_nearest(f):
    """Returns the int which is nearest to the argument float.  Draws
    are rounded away from zero. If the argument is None, returns
    None.
    """

    if f is None:
        return None
    if f >= 0:
        return int(math.floor(f + 0.5))
    else:
        return int(math.ceil(f - 0.5))

        


def read_antarc_line(line):
    """Convert a single line from the Antarctic/Australasian dataset
    files into a year and a 12-tuple of floats (the temperatures in
    Centigrade).  Missing values become None.
    """

    year = int(line[:4])
    line = line[4:]
    tuple = []
    if line[6] == '.' or line[7] == '-':
        # Some of the datasets are 12f8.1 with missing values as '       -'.
        for i in range(0,12):
            tuple.append(read_tenths(line[i*8:i*8+8]))
    else:
        # Others are xx12f7.1 or xxx12f7.1 with missing values as '       '.
        np = line.find('.')
        if np < 0:
            raise ValueError, "Non-data line encountered: '%s'" % line
        position = (np % 7) + 2
        for i in range(0,12):
            tuple.append(read_tenths(line[i*7+position:i*7+7+position]))
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


antarc_discard_re = re.compile(r'^$|^Get |^[12A-Z]$')
antarc_temperature_re = re.compile(r'^(.*) .* *temperature')


class USHCNReader(object):
    def __init__(self, path, station_path, record_discriminator=None):
        self.f = open_or_uncompress(path)
        self.record_discriminator = record_discriminator
        self.stations = read_USHCN_stations(station_path)
        self.us_only = {}
        for id12 in self.stations.itervalues():
            self.us_only[id12] = None

    def __iter__(self):
        return self._it()

    # TODO: This is reading data similar to other line based records,
    #       but the algorithm is somewhat different. See examples that
    #       use itertools.groupby.
    # NOTE: Sometimes, there are two sets of data for a station.
    #       Where the same year appears in both data sets, the second year's
    #       data is used. I do not know whether this is intentional or not.
    def _it(self):
        def fill_record():
            for year, temps in sorted(years_data.iteritems()):
                record.add_year_of_tenths(year, temps)
            return record

        record = None
        stations = self.stations
        years_data = {}
        for line in self.f:
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
                    # tenths of degree centigrade
                    temp = round_to_nearest((temp_fahrenheit - 320) * 5/9.0)
                    valid = True
                temps.append(temp)
            if valid: # some valid data found
                # We cannot add using 'add_year_of_tenths' because that breaks
                # things when years are duplicated in the data; (see not
                # above). The record is actually filled using ``fil_record``.
                #record.add_year_of_tenths(year, temps)
                years_data.setdefault(year, [])
                years_data[year] = temps

        if record is not None:
            yield fill_record()
        self.f.close()


def read_USHCN_stations(path):
    """Reads the USHCN station list and returns a dictionary
    mapping USHCN station ID to 12-digit station ID.
    """

    stations = {}
    for line in open(path):
        (USHCN_id, WMO_id, duplicate) = line.split()
        USHCN_id = int(USHCN_id)
        if WMO_id[0:3] != '425': # non-US country_code
            raise ValueError, "non-425 station found in ushcn.tbl: '%s'" % line
        if duplicate != '0':
            raise ValueError, "station in ushcn.tbl with non-zero duplicate: '%s'" % line
        stations[USHCN_id] = WMO_id + '0'
    # some USHCNv2 station IDs convert to USHCNv1 station IDs:
    for line in open('input/ushcnV2_cmb.tbl'):
        (v2_station,_,v1_station,_) = line.split()
        stations[int(v2_station)] = stations[int(v1_station)]
    return stations


def HohenpeissenbergReader(path):
    """reads the Hohenpeissenberg data from
    input/t_hohenpeissenberg_200306.txt_as_received_July17_2003
    which has a header line and then one line per year from 1781.
    We only want data from 1880 to 2002.
    """
    f = open(path)

    record = code.giss_data.StationRecord('617109620002')
    for line in f:
        if line[0] in '12':
            year = int(line[:4])
            if year < 1880 or year > 2002:
                continue
            data = line.split()
            temps = map(read_tenths, data[1:13])
            assert len(temps) == 12
            record.add_year_of_tenths(year, temps)
    yield record


def read_tenths(s):
    """Returns the integer nearest to the argument string times 10.
    If float conversion fails, returns MISSING.
    """

    try:
        f = float(s)
    except:
        return code.giss_data.MISSING
    return round_to_nearest(f * 10)



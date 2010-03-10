#!/usr/bin/env python
# $URL$
# $Rev$

"""Classes for GISTEMP data.

Primarily, the classes herein support temperature records, which
are composed of series of monthly averages. There two types of record,
derived from the `MonthlyTemperatureRecord` class:

`StationRecord`
    Stores a set of monthly averages associated with a particular monitoring
    `Station`.

`SubboxRecord`
    Stores combined monthly averages for StationRecords within a sub-box.

Both types of record can be grouped in collections, often (in the original
GISTEMP code) in files. Collections of records have associated metadata,
the `StationMetaData` and `SubboxMetaData` classes.

"""
__docformat__ = "restructuredtext"

import sys

import read_config

#: The base year for time series data. Data before this time is not
#: used in calculations.
BASE_YEAR = 1880

#: The value that is used to indicate a bad or missing data point.
MISSING = 9999.0

def invalid(v):
    return v == MISSING

def valid(v):
    return not invalid(v)


_stations = None
_v2_sources = None
_ghcn_last_year = None

def _load_v2_inv():
    """Load the v2.inv file.

    :Return:
        A dictionary of `Station` instances, keyed by the station's `uid`.

    """
    return dict((s.uid, s) for s in StationInventoryReader("input/v2.inv"))


def get_ghcn_last_year():
    """Get the latest year in the GHCN data.

    This simply reads the ``input/v2.mean`` file and extracts the year from
    each line.

    In the original GISTEMP code, this piece of information was cached in the
    file ``work/GHCN.last_year`` by step0. This alternative approach avoids the
    need for this file and perfectly quick enough.
    """
    global _ghcn_last_year
    if _ghcn_last_year is None:
        f = open("input/v2.mean")
        max_year = 0
        for l in f:
            if l[12:13] == '2': # first digit of year is a '2'
                max_year = max(int(l[12:16]), max_year)
        _ghcn_last_year = max_year
        f.close()

    return _ghcn_last_year


class StationInventoryReader(object):
    """Reader for files in the format of input/v2.inv.

    This can be used as a single-pass iterator, yielding `giss_data.Station`
    instances.

    For a list of metadata fields, see *fields*.

    The input file is in the same format as the GHCN V2 file v2.temperature.inv
    (in fact, it's the same file, but with records added for the Antarctic
    stations that GHCN doesn't have).  The best description of that file's
    format is the Fortran program:
    ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/v2/v2.read.inv.f

    Here are two typical lines, with a record diagram

    40371148001 ALMASIPPI,MA                    49.55  -98.20  274  287R   -9FLxxno-9x-9COOL FIELD/WOODSA1   0
    42572530000 CHICAGO/O'HARE, ILLINOIS        42.00  -87.90  205  197U 6216FLxxno-9A 1COOL CROPS      C3 125

    0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
    id---------xname--------------------------xlat---xlon----x1---2----34----5-6-7-8-910grveg-----------GU--11

       id                  40371148001          42572530000
       name                ALMASIPPI,MA         CHICAGO/O'HARE, ILLINOIS
       lat                 49.55                42.00
       lon                 -98.20               -87.90
    1  elevs               274                  205
    2  elevg               287                  197
    3  pop                 R                    U
    4  ipop                -9                   6216
    5  topo                FL                   FL
    6  stveg               xx                   xx
    7  stloc               no                   no
    8  iloc                -9                   -9
    9  airstn              x                    A
    10 itowndis            -9                   1
       grveg               COOL FIELD/WOODS     COOL CROPS
    G  GHCN_brightness     A                    C
    U  US_brightness       1                    3
    11 global_brightness   0                    125
    """
    fields = (
        (0,   11,  'uid',              str),
        (12,  42,  'name',             str),
        (43,  49,  'lat',              float),
        (50,  57,  'lon',              float),
        (58,  62,  'elevation',        int),
        (62,  67,  'ground_elevation', int),
        (67,  68,  'pop',              str),
        (68,  73,  'ipop',             int),
        (73,  75,  'topo',             str),
        (75,  77,  'stveg',            str),
        (77,  79,  'stloc',            str),
        (79,  81,  'iloc',             str),    #int
        (81,  82,  'airstn',           str),
        (82,  84,  'itowndis',         str),    #int
        (84,  100, 'grveg',            str),
        (100, 101, 'GHCN_brightness',  str),
        (101, 102, 'US_brightness',    str),
        (102, 106, 'global_brightness',int),
    )

    def __init__(self, path):
        """Constructor:

        :Param path:
            The path of the inventory file.

        """
        self.f = open(path)

    def __iter__(self):
        return self._it()

    def _it(self):
        for line in self.f:
            dd = dict((name, conv(line[a:b]))
                    for a, b, name, conv in self.fields)
            yield Station(**dd)
        self.f.close()


def stations():
    """Return a dictionary of all known stations.

    Note: This uses information from the ``input/v2.inv`` file. This file
    is loaded the first time this function is called. The time is reasonably
    short (around 0.3 second on a 2GHz celeron).

    :Return:
        A dictionary, keyed by a `Station.uid`. Each entry is a `Station`
        instance.

    """
    global _stations
    if _stations is None:
        _stations = _load_v2_inv()
    return _stations


def v2_sources():
    global _v2_sources
    if _v2_sources is None:
        _v2_sources = read_config.v2_get_sources()
    return _v2_sources


def clear_cache(func):
    """A decorator, for `TemperatureRecord` methods that change the data.

    Any method that changes the underlying data series in a `TemperatureRecord`
    must clear the cached values used for some properties. This decorator
    takes care of this chore.

    """
    def f(self, *args, **kwargs):
        self._good_count = None
        self._ann_anoms_good_count = None
        return func(self, *args, **kwargs)

    return f


class StationMetaData(object):
    """The metadata for a set of station records.

    :Ivar mo1:
        The number of months covered in the entire data set.
    :Ivar kq:
        The KQ quantity from the header record.
    :Ivar mavg:
        A code indicating the length of time each average value represents. The
        only supported by the CCC code is '6', which indicates that each entry
        is a monthly average. The effect of this having a different value is
        undefined.
    :Ivar monm:
        Maximum length of any time record.
    :Ivar monm4:
        This is the size of this record when written to a GISS Fortran
        unformatted file.

        TODO: This can probably be ditched and calculated as required
        in I/O code.
    :Ivar yrbeg:
        The year of first data.
    :Ivar missing_flag:
        The value used to indicate a missing value in the `series`. This is
        often referred to in other code as variously bad, BAD, XBAD.

        This should become unimportant over time in the CCC code, which should
        stick to always using the `MISSING` value.
    :Ivar precipitation_flag:
        Probably defines a special value that serves a similar purppose to
        the `missing_flag`. This does not seem to be used by any CCC code.
    :Ivar mlast:
        TODO
    :Ivar title:
        A title for this set of station records.
    """
    def __init__(self, **k):
        self.__dict__ = k

    def __repr__(self):
        return 'StationMetadata(%r)' % self.__dict__


class Station(object):
    """A monitoring station's information.

    This holds the information about a single monitoring station. Not all the
    fields are used by the CCC code and (currently) the meaning of some of
    these is unknown.

    :Ivar lat, lon:
        The location of the station as floats, representing degrees.
    :Ivar elevation:
        The station's height in metres.
    :Ivar uid:
        The unique ID of the station. This is held as an 11 digit string.
    :Ivar name:
        The station's name.
    :Ivar GHCN_brightness:
        TODO: Some indication of brightness.
    :Ivar US_brightness:
        Brighness indication for US stations. A value of '1' or ' ' is taken
        to indicate that the station is in a rural area.
    :Ivar airstn:
        Unknown.
    :Ivar ground_elevation, elevation:
        The ground elevation at the station's location and
        the elevation of the station above the ground level. Both in metres.
    :Ivar grveg:
        An indication of the type of ground vegetation. For example,
        'TROPICAL DRY FOR'.
    :Ivar stveg:
        Unknown.
    :Ivar iloc:
        Unknown.
    :Ivar ipop:
        Unknown.
    :Ivar itowndis:
        Unknown.
    :Ivar pop:
        Unknown.
    :Ivar stloc:
        Unknown.
    :Ivar topo:
        Unknown.
    :Ivar global_brightness:
        A global brightness index, range 0-186 (at least)

    """
    def __init__(self, name=None, lat=None, lon=None, uid=None,
            ground_elevation=None,
            elevation=None, pop=None, ipop=None, topo=None, stveg=None,
            stloc=None, iloc=None, airstn=None, itowndis=None,
            grveg=None, GHCN_brightness=None, US_brightness=None,
            global_brightness=None):
        self.name = name
        self.lat = lat
        self.lon = lon
        self.uid = uid
        self.ground_elevation = ground_elevation
        self.elevation = elevation
        self.pop = pop
        self.ipop = ipop
        self.topo = topo
        self.stveg = stveg
        self.stloc = stloc
        self.iloc = iloc
        self.airstn = airstn
        self.itowndis = itowndis
        self.grveg = grveg
        self.GHCN_brightness = GHCN_brightness
        self.US_brightness = US_brightness
        self.global_brightness = global_brightness


# TODO: Needs some review. Among things to think about:
#
# 1. Might it be seen as too complicated? It is complicated for a reason; to
#    make the code that manipulates temperature series more readable.
# 2. Should we use properties or convert the properties to methods?
# 3. Some of the names are open to improvement.
class MonthlyTemperatureRecord(object):
    """Base class for monthly temperature records.

    This is the base class for both the `StationRecord` and `SubboxRecord`
    classes. It contains a series of average monthly temperatures, which are
    accessible via the `series` property. The series propery always provides
    an array of floating point values in celsius. This property should
    **always** be treated as read-only; the effect of modifying elements is
    undefined.

    The series coveres the months from `first_month` to `last_month` month
    inclusive. Months are counted from a non-existant year zero. So January,
    1 AD has a month number of 13, February is 14, etc.

    Within series, some leading months and some trailing months may be set to
    the `MISSING` value. The `good_start_idx` and `good_end_idx` members define
    the Python range within the series that excludes these missing value.

    The GISTEMP/CCC code only uses data that starts from `BASE_YEAR` (1880).
    Some code works on data series that start from this base year. So it is
    convenient to be able to work in terms of years and months relative to this
    base year. There are a number of properties with names that start with
    `rel_` that provide values using this alternative reference.

    Note that most of the series metadata is provided by properties, which
    are effectively read-only. All the instance variables should also be
    treated as read-only and you should only set values in the data series
    using the provided methods.

    """
    def __init__(self):
        self._first_month = sys.maxint
        self._good_start_idx = sys.maxint
        self._good_end_idx = 0
        self._series = []
        self._good_count = None
        self._ann_anoms_good_count = None

    def is_empty(self):
        """Test whether the record contains data."""
        return len(self._series) == 0

    @property
    def n(self):
        """The length of the series."""
        return len(self._series)

    @property
    def good_start_idx(self):
        """Index of the first good value in the `series`.

        It is always true that ``series[good_start_idx:good_end_idx]`` will
        either be empty or start and end with a valid value.

        """
        return self._good_start_idx

    @property
    def good_end_idx(self):
        """Index of the entry after the last good value in the `series`.

        It is always true that ``series[good_start_idx:good_end_idx]`` will
        either be empty or start and end with a valid value.

        """
        return self._good_end_idx

    @property
    def first_month(self):
        """The number of the last months in the data series.

        This number is counted from January in a non-existant year zero. The
        property `last_month` provides the other end of the inclusive range of
        months held in the `series`.

        The `series` contains `last_month` - `first_month` + 1 entries.

        """
        return self._first_month

    @property
    def last_month(self):
        """The number of the last months in the data series.

        The `series` contains ``last_month`` - `first_month` + 1 entries.
        See `first_month` for details of how months are counted.

        """
        return (self.first_month + len(self._series) - 1)

    @property
    def first_year(self):
        """The year of the first value in the series."""
        return (self.first_month - 1) // 12

    @property
    def last_year(self):
        """The year of the last value in the series."""
        return (self.last_month - 1) // 12

    @property
    def first_good_year(self):
        """The year of the first good value in the series."""
        return (self.first_good_month - 1) // 12

    @property
    def last_good_year(self):
        """The year of the last good value in the series."""
        return (self.last_good_month - 1) // 12

    @property
    def first_good_month(self):
        """The number of the first good month in the series.

        See `first_month` for details of how months are counted.

        """
        return self.first_month + self.good_start_idx

    @property
    def last_good_month(self):
        """The number of the last good month in the series.

        See `first_month` for details of how months are counted.

        """
        return self.first_month + self.good_end_idx - 1

    @property
    def rel_first_year(self):
        """The `first_year` relative to `BASE_YEAR`."""
        return self.first_year - BASE_YEAR

    @property
    def rel_last_year(self):
        """The `last_year` relative to `BASE_YEAR`."""
        return self.last_year - BASE_YEAR

    @property
    def rel_first_good_year(self):
        """The `first_good_year` relative to `BASE_YEAR`."""
        return self.first_good_year - BASE_YEAR

    @property
    def rel_last_good_year(self):
        """The `last_good_year` relative to `BASE_YEAR`."""
        return self.last_good_year - BASE_YEAR

    @property
    def rel_first_month(self):
        """The `first_month` relative to `BASE_YEAR`."""
        return self.first_month - BASE_YEAR * 12

    @property
    def rel_last_month(self):
        """The `last_month` relative to `BASE_YEAR`."""
        return self.last_month - BASE_YEAR * 12

    @property
    def rel_first_good_month(self):
        """The `first_good_month` relative to `BASE_YEAR`."""
        return self.first_good_month - BASE_YEAR * 12

    @property
    def rel_last_good_month(self):
        """The `last_good_month` relative to `BASE_YEAR`."""
        return self.last_good_month - BASE_YEAR * 12

    @property
    def good_count(self):
        """The number of good values in the data."""
        if self._good_count is None:
            bad_count = 0
            for v in self._series:
                bad_count += invalid(v)
            self._good_count = len(self._series) - bad_count
        return self._good_count

    def strip_invalid(self):
        """Strip leading and trailing invalid values.

        Adjusts the record so that the series starts and ends with a good (not
        `MISSING`) value. If there are no good values, the series will be
        emptied.

        """
        self._first_month = self.first_good_month
        self._series[:] = self._series[self._good_start_idx:self._good_end_idx]
        self._good_start_idx = 0
        self._good_end_idx = len(self._series)
    strip_invalid = clear_cache(strip_invalid)

    def get_monthly_valid_counts(self):
        """Get number of good values for each month.

        :Return:
            A list of 12 entries. Entry zero is the number of good entries
            for January, entry 1 for february, etc.

        """
        monthly_valid = [0] * 12
        for i, v in enumerate(self._series):
            monthly_valid[(self.first_month + i - 1) % 12] += valid(v)
        return monthly_valid

    def _set_series(self, first_month, series):
        self._first_month = first_month
        self._good_start_idx = sys.maxint
        self._good_end_idx = 0
        self._series = []
        for v in series:
            if invalid(v):
                self._series.append(MISSING)
            else:
                self._good_start_idx = min(self._good_start_idx,
                        len(self._series))
                self._series.append(v)
                self._good_end_idx = max(self._good_end_idx, len(self._series))
    _set_series = clear_cache(_set_series)

    def add_year(self, year, data):
        if self.first_month != sys.maxint:
            # We have data already, so we may need to pad with missing months
            # Note: This assumes the series is a whole number of years.
            gap = year - self.last_year - 1
            if gap > 0:
                self._series.extend([MISSING] * gap * 12)
        start_month = year * 12 + 1
        self._first_month = min(self.first_month, start_month)
        for v in data:
            if invalid(v):
                self._series.append(MISSING)
            else:
                self._good_start_idx = min(self._good_start_idx,
                        len(self._series))
                self._series.append(v)
                self._good_end_idx = max(self._good_end_idx, len(self._series))
    add_year = clear_cache(add_year)


class StationRecord(MonthlyTemperatureRecord):
    """An average monthly temperature record associated with a `Station`.

    There can be multiple temperature series for a single `Station`. The
    `station` property provides the associated `Station` instance.

    :Ivar uid:
        An integer that acts as a unique ID for the time series. This
        is generated by taking the last 8 digits of the `Station` uid, and
        adding an additional digit. The last digit distinguished this series
        from other series from the same station.
    :Ivar source:
        The source of the data, which defaults to 'UNKNOWN'.

    """
    def __init__(self, uid, **kwargs):
        super(StationRecord, self).__init__()
        self.uid = uid
        self.source = v2_sources().get(uid, "UNKNOWN")
        self.ann_anoms = []

    def __str__(self):
        return "StationRecord(uid=%r)" % self.uid

    @property
    def series(self):
        """The series of values in celsius."""
        return self._series

    @property
    def station(self):
        """The corresponding `Station` instance."""
        st = stations().get(self.station_uid)
        if st is None:
            print "BUM!", self.uid, self.station_uid
        return st

    @property
    def station_uid(self):
        """The unique ID of the corresponding station."""
        try:
            return self.uid[:-1]
        except:
            print ">>>", self, self.uid
            raise

    @property
    def short_id(self):
        """The shortened form of the record's uinique ID, as an integer.

        This is the `uid` with the country code removed, converted to an
        integer.

        """
        return int(self.uid[3:])

    @property
    def country_code(self):
        """The country code as an integer.

        This is the first three characters of the `uid`, converted to an
        integer.

        """
        return int(self.uid[:3])

    @property
    def discriminator(self):
        """The discriminator code.

        This distingushes different records associated with same station.
        This is the last digit of the record's `uid` converted to an
        integer.

        """
        return int(self.uid[-1])

    def has_data_for_year(self, year):
        for t in self.get_a_year(year):
            if t != MISSING:
                return True

    def get_a_month(self, month):
        """Get the value for a single month."""
        idx = month - self.first_month
        if idx < 0:
            return MISSING
        try:
            return self.series[month - self.first_month]
        except IndexError:
            return MISSING

    def get_a_year(self, year):
        """Get the time series data for a year."""
        start_month = year * 12 + 1
        return [self.get_a_month(m) for m in range(start_month, start_month + 12)]

    def get_set_of_years(self, first_year, last_year):
        """Get a set of year records.

        :Return:
            A list of lists, where each sub-list contains 12 temperature values
            for a given year. This works for any range of years, missing years
            are filled with the MISSING value.

        """
        return [self.get_a_year(y) for y in range(first_year, last_year + 1)]

    def set_series(self, first_month, series):
        self._set_series(first_month, series)

    def set_ann_anoms(self, ann_anoms):
        self.ann_anoms[:] = ann_anoms

    @property
    def ann_anoms_good_count(self):
        """Number of good values in the annual anomalies"""
        if self._ann_anoms_good_count is None:
            bad = 0
            for v in self.ann_anoms:
                bad += invalid(v)
            self._ann_anoms_good_count = len(self.ann_anoms) - bad
        return self._ann_anoms_good_count

    def copy(self):
        r = StationRecord(self.uid)
        r.set_series(self.first_month, self.series)
        return r

    def report_str(self):
        """Format neatly for reporting purposes."""

        s = "StationRecord: %s - %s\n" % (self.station.name, self.uid)
        s += "   All data  : %d:%02d - %d:%02d  [%3d:%02d - %3d:%02d]" % (
                self.first_year, (self.first_month - 1) % 12 + 1,
                self.last_year, (self.last_month  - 1) % 12 + 1,
                self.rel_first_year, (self.first_month - 1) % 12 + 1,
                self.rel_last_year, (self.last_month  - 1) % 12 + 1,
                )
        s += "  %5d - %5d  (%4d - %4d)\n" % (
                self.first_month, self.last_month,
                self.rel_first_month, self.rel_last_month,
                )
        s += "   Good data : %d:%02d - %d:%02d  [%3d:%02d - %3d:%02d]  %d-%d" % (
                self.first_good_year, (self.first_good_month - 1) % 12 + 1,
                self.last_good_year, (self.last_good_month - 1) % 12 + 1,
                self.rel_first_good_year, (self.first_good_month - 1) % 12 + 1,
                self.rel_last_good_year, (self.last_good_month - 1) % 12 + 1,
                self.good_start_idx, self.good_end_idx,
                )
        s += "  %5d - %5d  (%4d - %4d)\n" % (
                self.first_good_month, self.last_good_month,
                self.rel_first_good_month, self.rel_last_good_month,
                )
        for n in ("series", "anomalies"):
            v = getattr(self, n, None)
            if v is None:
                continue
            s += "   %-10s: [%-4d] = %s...\n" % (n, len(v),
                    str([str(x)[:6] for x in v[:6]])[:-1])
            s += "   %-10s:            ...%s\n" % ("",
                    str([str(x)[:6] for x in v[-6:]])[1:])

        return s[:-1]


class SubboxMetaData(object):
    """The metadata for a set of sub-box records.

    :Ivar mo1:
       TBD
    :Ivar kq:
       TBD
    :Ivar mavg:
       TBD
    :Ivar monm:
       TBD
    :Ivar monm4:
       TBD
    :Ivar yrbeg:
       TBD
    :Ivar missing_flag:
       TBD
    :Ivar precipitation_flag:
       TBD
    :Ivar title:
       TBD
    """
    def __init__(self, mo1, kq, mavg, monm, monm4, yrbeg,
            missing_flag, precipitation_flag, title):
        self.mo1 = mo1
        self.kq = kq
        self.mavg = mavg
        self.monm = monm
        self.monm4 = monm4
        self.yrbeg = yrbeg
        self.missing_flag = missing_flag
        self.precipitation_flag = precipitation_flag
        self.title = title

    def __repr__(self):
        return '<StationMetadata %r>' % self.__dict__

    def copy(self):
        return SubboxMetaData(self.mo1, self.kq, self.mavg, self.monm,
                self.monm4, self.yrbeg, self.missing_flag,
                self.precipitation_flag, self.title)


class SubboxRecord(MonthlyTemperatureRecord):
    """A sub-box record.

    This can hold, for example, a record of data as stored within the
    ``input/SBBX.HadR2`` file.

    :Ivar lat_S, lat_N, lon_W, lon_E:
        Coordinates describing the box's area.
    :Ivar stations:
        The number of stations that contributed to this sub-box.
    :Ivar station_months:
        The number of months that contributed to this sub-box.
    :Ivar d:
        Unknown.

    """
    def __init__(self, lat_S, lat_N, lon_W, lon_E,
            stations, station_months, d, series, **_ignored):
        super(SubboxRecord, self).__init__()
        self.lat_S = lat_S
        self.lat_N = lat_N
        self.lon_W = lon_W
        self.lon_E = lon_E
        self.stations = stations
        self.station_months = station_months
        self.d = d
        self.set_series(series)

    @property
    def series(self):
        """The series of values in celsius."""
        return self._series

    def set_series(self, series):
        self._set_series(BASE_YEAR * 12, series)

    def pad_with_missing(self, n):
        while self.n < n:
            self._series.append(MISSING)
    pad_with_missing = clear_cache(pad_with_missing)

    def set_value(self, idx, value):
        while idx >= len(self.series):
            self._series.append(MISSING)
        self._series[idx] = value
    set_value = clear_cache(set_value)

    def trim(self):
        self.station_months = len(self._series) - self._series.count(
                MISSING)

    def __repr__(self):
        return ('<Subbox (%+06.2f,%+06.2f) (%+07.2f,%+07.2f): %d>' %
                (self.lat_S, self.lat_N, self.lon_W, self.lon_E, self.n))
    def report_str(self):
        """Format neatly for reporting purposes."""

        s = "SubboxRecord:\n"
        s += "   All data  : %d:%02d - %d:%02d  [%3d:%02d - %3d:%02d]" % (
                self.first_year, (self.first_month - 1) % 12 + 1,
                self.last_year, (self.last_month  - 1) % 12 + 1,
                self.rel_first_year, (self.first_month - 1) % 12 + 1,
                self.rel_last_year, (self.last_month  - 1) % 12 + 1,
                )
        s += "  %5d - %5d  (%4d - %4d)\n" % (
                self.first_month, self.last_month,
                self.rel_first_month, self.rel_last_month,
                )
        s += "   Good data : %d:%02d - %d:%02d  [%3d:%02d - %3d:%02d]  %d-%d" % (
                self.first_good_year, (self.first_good_month - 1) % 12 + 1,
                self.last_good_year, (self.last_good_month - 1) % 12 + 1,
                self.rel_first_good_year, (self.first_good_month - 1) % 12 + 1,
                self.rel_last_good_year, (self.last_good_month - 1) % 12 + 1,
                self.good_start_idx, self.good_end_idx,
                )
        s += "  %5d - %5d  (%4d - %4d)\n" % (
                self.first_good_month, self.last_good_month,
                self.rel_first_good_month, self.rel_last_good_month,
                )
        for n in ("series", "anomalies"):
            v = getattr(self, n, None)
            if v is None:
                continue
            s += "   %-10s: [%-4d] = %s...\n" % (n, len(v),
                    str([str(x)[:6] for x in v[:6]])[:-1])
            s += "   %-10s:            ...%s\n" % ("",
                    str([str(x)[:6] for x in v[-6:]])[1:])

        return s[:-1]

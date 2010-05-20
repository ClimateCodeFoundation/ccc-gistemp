#!/usr/bin/env python
# $URL$
# $Rev$

"""Classes for GISTEMP data.

Primarily, the classes herein support monthly temperature series.
Typically these are either for stations (station records) or subboxes
(subbox series).  In either case the same class is used, `Series`,
differing in what keyword arguments are supplied.

station records
    Stores a set of monthly averages associated with a particular monitoring
    `Station`.

subbox series
    Stores monthly averages for a subbox, typically synthesized by
    combining several station records together.

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


class Station(object):
    """A monitoring station's information.

    This holds the information about a single monitoring station. Not
    all the attributes are used by the CCC code.  For a list of
    attributes and documentation, see the stations() function.
    """
    def __init__(self, **values):
        self.__dict__.update(values)


def stations():
    """Return a dictionary of all known Stations, keyed by
    Station.uid.  The first time this is called, it reads the station
    inventory file input/v2.inv.

    The input file is in the same format as the GHCN V2 file v2.temperature.inv
    (in fact, it's the same file, but with records added for the Antarctic
    stations that GHCN doesn't have).   Descriptions of that file's
    format can be found in the Fortran programs:
    ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/v2/v2.read.inv.f
    ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/v2/v2.read.data.f

    Here are two typical lines, with a record diagram

    id---------xname--------------------------xlat---xlon----x1---2----34----5-6-7-8-910grveg-----------GU--11
    0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
    40371148001 ALMASIPPI,MA                    49.55  -98.20  274  287R   -9FLxxno-9x-9COOL FIELD/WOODSA1   0
    42572530000 CHICAGO/O'HARE, ILLINOIS        42.00  -87.90  205  197U 6216FLxxno-9A 1COOL CROPS      C3 125

       uid                 40371148001          42572530000
          The unique ID of the station. This is held as an 11 digit string.
       name                ALMASIPPI,MA         CHICAGO/O'HARE, ILLINOIS
        The station's name.
       lat                 49.55                42.00
        The latitude, in degrees to two decimal places.
       lon                 -98.20               -87.90
        The longitude, in degrees to two decimal places.
    1  elevs               274                  205
        The station elevation? in metres?
    2  elevg               287                  197
        The ground elevation? in metres?
    3  pop                 R                    U
        'R' for rural,  'S' for semi-urban, 'U' for urban
    4  ipop                -9                   6216
    5  topo                FL                   FL
        The topography
    6  stveg               xx                   xx
    7  stloc               no                   no
    8  iloc                -9                   -9
    9  airstn              x                    A
    10 itowndis            -9                   1
       grveg               COOL FIELD/WOODS     COOL CROPS
        An indication of the type of ground vegetation. For example,
        'TROPICAL DRY FOR'.
    G  GHCN_brightness     A                    C
    U  US_brightness       1                    3
        Brighness indication for US stations, or ' ' for non-US
        stations.  '1' is dark, '3' is bright.
    11 global_brightness   0                    125
        A global brightness index, range 0-186 (at least)
    """

    global _stations
    if _stations is None:
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
            (79,  81,  'iloc',             str),    #int?
            (81,  82,  'airstn',           str),
            (82,  84,  'itowndis',         str),    #int?
            (84,  100, 'grveg',            str),
            (100, 101, 'GHCN_brightness',  str),
            (101, 102, 'US_brightness',    str),
            (102, 106, 'global_brightness',int),
        )

        _stations = {}
        for line in open("input/v2.inv"):
            dd = dict((name, conv(line[a:b]))
                      for a, b, name, conv in fields)
            _stations[dd['uid']] = Station(**dd)

    return _stations



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


# TODO: Needs some review. Among things to think about:
#
# 1. Might it be seen as too complicated? It is complicated for a reason; to
#    make the code that manipulates temperature series more readable.
# 2. Should we use properties or convert the properties to methods?
# 3. Some of the names are open to improvement.
class Series(object):
    """Monthly temperature Series.

    Instances contain a series of average monthly temperatures, which are
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

    There are no subclasses of this class.  Some instances represent
    station records, other instances represent subbox series.

    For station records there can be multiple series for a single `Station`.
    The `station` property provides the associated `Station` instance.
    For a given station the different series are called "duplicates" in
    GHCN terminology; they have a 12-digit uid that is made up of an
    11-digit station identifier and a single extra digit to distinguish
    each of the station's series.

    Generally a station record will have its uid supplied as a keyword
    argument to the constructor (accessing the `station` property relies
    on this):

    :Ivar uid:
        An integer that acts as a unique ID for the time series. This
        is generally a 12-digit identifier taken from the GHCN file; the
        first 11 digits comprise an identifier for the station.
	The last digit distinguishes this series from other series
	from the same station.

    When used to hold a series for a subbox, for example a record of data
    as stored in the ``input/SBBX.HadR2`` file, then the following
    keyword arguments are traditionally supplied to the constructor:

    :Ivar lat_S, lat_N, lon_W, lon_E:
        Coordinates describing the box's area.
    :Ivar stations:
        The number of stations that contributed to this sub-box.
    :Ivar station_months:
        The number of months that contributed to this sub-box.
    :Ivar d:
        Characteristic distance to station closest to centre.

    """
    def __init__(self, **k):
        self._first_month = sys.maxint
        self._good_start_idx = sys.maxint
        self._good_end_idx = 0
        self._series = []
        self._good_count = None
        self._ann_anoms_good_count = None
        self.ann_anoms = []
        series = None
        if 'series' in k:
            series = k['series']
            del k['series']
            self.set_series(BASE_YEAR*12+1, series)
        self.__dict__.update(k)

        if hasattr(self, 'uid'):
            # Generally applies to station records
            self.source = v2_sources().get(self.uid, "UNKNOWN")
        elif hasattr(self, 'box'):
            # Generally applies to subbox series.
            # Synthesize a uid attribute based on the box's centre.
            import eqarea
            lat,lon = eqarea.centre(self.box)
            self.uid = "%+05.1f%+06.1fC" % (lat,lon)

    def __repr__(self):
        # A bit ugly, because it tries to do something sensible for both
        # station records and subbox series.
        if hasattr(self, 'box'):
            return ('Series(box=(%+06.2f,%+06.2f,%+07.2f,%+07.2f))' %
              self.box)
        else:
            # Assume it is a station record with a uid.
            return "Series(uid=%r)" % self.uid

    @property
    def series(self):
        """The series of values (conventionally in degrees Celsius)."""
        return self._series

    def __len__(self):
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
        """The number of the first month in the data series.

        This number is counted from January (being 1) in a non-existant
	year zero. The property `last_month` provides the other end
	of the inclusive range of months held in the `series`.

        The `series` contains `last_month` - `first_month` + 1 entries.

        """
        return self._first_month

    @property
    def last_month(self):
        """The number of the last month in the data series.

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

    def get_monthly_valid_counts(self):
        """Get number of good values for each month.

        :Return:
            A list of 12 entries. Entry zero is the number of good entries
            for January, entry 1 for February, etc.

        """
        monthly_valid = [0] * 12
        for i, v in enumerate(self._series):
            monthly_valid[(self.first_month + i - 1) % 12] += valid(v)
        return monthly_valid

    # Year's worth of missing data
    missing_year = [MISSING]*12

    def has_data_for_year(self, year):
        return self.get_a_year(year) != self.missing_year

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
        start_index = start_month - self.first_month
        data = self.series[start_index:start_index+12]
        if len(data) == 12:
            return data
        # Do it the slow way:
        return [self.get_a_month(m)
                for m in range(start_month, start_month + 12)]

    def get_set_of_years(self, first_year, last_year):
        """Get a set of year records.

        :Return:
            A list of lists, where each sub-list contains 12 temperature values
            for a given year. This works for any range of years, missing years
            are filled with the MISSING value.

        """
        return [self.get_a_year(y) for y in range(first_year, last_year + 1)]

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

    def pad_with_missing(self, n):
        while len(self) < n:
            self._series.append(MISSING)
    pad_with_missing = clear_cache(pad_with_missing)

    def trim(self):
        self.station_months = len(self._series) - self._series.count(
                MISSING)

    @property
    def station(self):
        """The corresponding `Station` instance.  Only works for station
        records."""
        st = stations().get(self.station_uid)
        if st is None:
            print "BUM!", self.uid, self.station_uid
        return st

    @property
    def station_uid(self):
        """The unique ID of the corresponding station."""
        return self.uid[:11]

    # Mutators below here

    def set_series(self, first_month, series):
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

    def add_year(self, year, data):
        if self.first_month == sys.maxint:
            self._first_month = year * 12 + 1
        else:
            # We have data already, so we may need to pad with missing months
            # Note: This assumes the series is a whole number of years.
            gap = year - self.last_year - 1
            if gap > 0:
                self._series.extend([MISSING] * gap * 12)
        assert self.first_month % 12 == 1
        if year < self.first_year:
            # Ignore years before the first year.  Previously this case
            # was extremely buggy.
            print self.uid, year, self.first_year
            return
        assert year == self.last_year + 1
         
        for v in data:
            if invalid(v):
                self._series.append(MISSING)
            else:
                self._good_start_idx = min(self._good_start_idx,
                        len(self._series))
                self._series.append(v)
                self._good_end_idx = max(self._good_end_idx, len(self._series))
    add_year = clear_cache(add_year)

    def set_value(self, idx, value):
        while idx >= len(self.series):
            self._series.append(MISSING)
        self._series[idx] = value
    set_value = clear_cache(set_value)


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
        return 'SubboxMetadata(%r)' % self.__dict__


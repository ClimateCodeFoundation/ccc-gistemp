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
import copy
import math

#: The base year for time series data. Data before this time is not
#: used in calculations.
BASE_YEAR = 1880

#: Integer code used to indicate missing data. This is units of 0.1
#: celcius.
MISSING = 9999

#: The floating point version of `MISSING`.
XMISSING = float(MISSING)

_stations = None

def _load_v2_inv():
    """Load the v2.inv file.

    :Return:
        A dictionary of `Station` instances, keyed by the station's
        `uid`.

    """
    return dict((s.uid, s) for s in StationInventoryReader("input/v2.inv"))


class StationInventoryReader(object):
    """Reader for files in the format of input/v2.inv.
    
    This can be used as a single-pass iterator, yielding `giss_data.Station`
    instances.

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
        (102, 106, 'idontknow',        str),
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


# TODO: This load-on demand approach is less than perfect because the loading
#       takes seconds.
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


def clear_cache(func):
    """A decorator, for `TemperatureRecord` methods that change the data.

    Any method that changes the underlying data series in a `TemperatureRecord`
    must clear the cached values used for some properties. This decorator
    takes care of this chore.

    """
    def f(self, *args, **kwargs):
        self._tenths = None
        self._celcius = None
        self._ngood = None
        self._ngood_ann_anoms = None
        return func(self, *args, **kwargs)

    return f


class StationMetaData(object):
    """The metadata for a set of station records.

    TODO: The descriptions of the header entries are work in progress.
    
    :Ivar mo1:
       TBD
    :Ivar kq:
        The KQ quantity from the header record.
    :Ivar mavg:
        A code indicating the lenght of time each average value represents. The
        only supported value is '6', which indicates that each entry is a
        monthly average. The effect of this having a different value is
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
        TBD
    :Ivar title:
        A title for this set of station records.
    """
    def __init__(self, mo1, kq, mavg, monm, monm4, yrbeg,
            missing_flag, precipitation_flag, mlast, title,):
        self.mo1 = mo1
        self.kq = kq
        self.mavg = mavg
        self.monm = monm
        self.monm4 = monm4
        self.yrbeg = yrbeg
        self.missing_flag = missing_flag
        self.precipitation_flag = precipitation_flag
        self.mlast = mlast
        self.title = title

    def __repr__(self):
        return '<StationMetadata %r>' % self.__dict__

    # TODO: This may be slightly dubious, but the code that uses looks fairly
    #       clean. It might be better to either have a factory function to do
    #       the conversion.
    def make_sub_box_meta(self):
        """Create a SubboxMetaData instance with matching values.

        """
        return SubboxMetaData(self.mo1, self.kq, self.mavg, self.monm,
                self.monm + 7, self.yrbeg, self.missing_flag,
                self.precipitation_flag, self.title)

    def pretty(self):
        """Format prettily, only really for debugging."""
        s = "StationMetaData:\n"
        for n, v in sorted(self.__dict__.iteritems()):
            if not n.startswith("_"):
                s += "  %-20s: %r\n" % (n, v)
        return s[:-1]
           

class Station(object):
    """A monitoring station's information.
    
    This holds the information about a single monitoring station.

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
        TODO
    :Ivar ground_elevation, elevation:
        The ground elevation at the station's location and
        the elevation of the station above the ground level. Both in metres.
    :Ivar grveg:
        An indication of the type of ground vegetation. For example,
        'TROPICAL DRY FOR'.
    :Ivar stveg:
        TODO: Something about the station's vegitation?
    :Ivar iloc:
        TODO
    :Ivar ipop:
        TODO
    :Ivar itowndis:
        TODO
    :Ivar pop:
        TODO
    :Ivar stloc:
        TODO
    :Ivar topo:
        TODO
    :Ivar idontknow:
        The last 3 digits in each line of v2.inv.

    """
    def __init__(self, name=None, lat=None, lon=None, uid=None,
            ground_elevation=None,
            elevation=None, pop=None, ipop=None, topo=None, stveg=None,
            stloc=None, iloc=None, airstn=None, itowndis=None,
            grveg=None, GHCN_brightness=None, US_brightness=None,
            idontknow=None):
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
        self.idontknow = idontknow

    @property
    def lat_fixed_1(self):
        """The latitude rounded to 1 decimal place."""
        return 0.1 * math.floor(self.lat * 10 + 0.5)

    @property
    def lon_fixed_1(self):
        """The longitude rounded to 1 decimal place."""
        return 0.1 * math.floor(self.lon * 10 + 0.5)

    @property
    def lat_as_tenths(self):
        """The latitude as a integer number of 0.1 degrees."""
        return int(math.floor(self.lat * 10 + 0.5))

    @property
    def lon_as_tenths(self):
        """The longitude as a integer number of 0.1 degrees."""
        return int(math.floor(self.lon * 10 + 0.5))

    def pretty(self):
        """Format prettily, only really for debugging."""
        s = "Station:\n"
        for n, v in sorted(self.__dict__.iteritems()):
            if n.startswith("_"):
                continue
        return s[:-1]

           
# TODO: Needs some review. Among things to think about:
#
# 1. Might it be seen as too complicated? It is complicated for a reason; to
#    make the code that manipulates temperature series more readable.
# 2. Should we use properties or convert the properties to methods?
# 3. Should all instance variables become properties (or methods)?
# 4. Some of the names are open to improvement.
class MonthlyTemperatureRecord(object):
    """Base class for monthly temperature records.

    This is the base class for both the `StationRecord` and `SubboxRecord`
    classes. It contains a series of average monthly temperatures, which are
    accessible via the `series` property. The series propery always provides
    an array of floating point values in celcius. This property should
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

    :Ivar good_start_idx, good_end_idx:
        The range of values in `series`. So ``series[good_start_idx:good_end_idx]``
        will either be empty or start and end with a valid value.
    :Ivar series:
        The temperature series for this station. The values are stored
        in celsius.
        (TODO currently sometimes tenths of a degree, but
        that needs to be changed.)
    :Ivar first_month:
        The number of the first month in the data series, counting from a
        January in a non-existant year zero. The property `last_month` provides
        the other end of the inclusive range of months held in the `series`.

    """
    def __init__(self):
        self.first_month = sys.maxint
        self.good_start_idx = sys.maxint
        self.good_end_idx = 0
        self._series = []
        self._celcius = None
        self._tenths = None
        self._ngood = None
        self._ngood_ann_anoms = None

    @classmethod
    def valid(cls, v):
        return not cls.invalid(v)

    @property
    def n(self):
        """The length of the series."""
        return len(self._series)

    @property
    def last_month(self):
        """The number of the last months in the data series.

        The `series` contains ``last_month`` - `first_month` + 1 entries.

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
        """TODO"""
        return self.first_month + self.good_start_idx

    @property
    def last_good_month(self):
        """TODO"""
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
    def ngood(self):
        """The number of good values in the data."""
        # TODO: Rename or remove or something. Also handle magic No.
        if self._ngood is None:
            bad = 0
            for v in self._series:
                bad += self.invalid(v)
            self._ngood = len(self._series) - bad
        return self._ngood

    @clear_cache
    def strip_invalid(self):
        """Strip leading and trailing invalid values."""
        self.first_month = self.first_good_month
        self._series[:] = self._series[self.good_start_idx:self.good_end_idx]
        self.good_start_idx = 0
        self.good_end_idx = len(self._series)

    def get_monthly_valid_counts(self):
        monthly_valid = [0] * 12
        for i, v in enumerate(self._series):
            monthly_valid[(self.first_month + i - 1) % 12] += self.valid(v)
        return monthly_valid

    @clear_cache
    def _set_series(self, first_month, series, missing, convert=lambda x:x):
        self.first_month = first_month
        self.good_start_idx = sys.maxint
        self.good_end_idx = 0
        self._series[:] = []
        for in_value in series:
            v = convert(in_value)
            if self.invalid(v):
                self._series.append(missing)
            else:
                self.good_start_idx = min(self.good_start_idx, len(self._series))
                self._series.append(v)
                self.good_end_idx = max(self.good_end_idx, len(self._series))

    @clear_cache
    def _add_year_of_data(self, year, data, missing, convert=lambda x:x):
        if self.first_month != sys.maxint:
            # We have data already, so we may need to pad with missing months
            # TODO: This is not fully correct because it assumes the series is
            #       already a whole number of years.
            gap = year - self.last_year - 1
            if gap > 0:
                self._series.extend([missing] * gap * 12)
        start_month = year * 12 + 1
        self.first_month = min(self.first_month, start_month)
        for m, in_value in enumerate(data):
            v = convert(in_value)
            if self.invalid(v):
                self._series.append(missing)
            else:
                self.good_start_idx = min(self.good_start_idx, len(self._series))
                self._series.append(v)
                self.good_end_idx = max(self.good_end_idx, len(self._series))


class StationRecord(MonthlyTemperatureRecord):
    """An average monthly temperature record associated with a `Station`.
    
    There can be multiple temperature series for a single `Station`. The
    `station` property provides the associated `Station` instance.

    :Ivar country_code:
        The first three digits of the corresponding `Station` uid.
    :Ivar uid:
        An integer that acts as a unique ID for the time series. This
        is generated by taking the last 8 digits of the `Station` uid, and
        adding an additional digit. The last digit distinguished this series
        from other series from the same station.
    :Ivar years:
        The number of years for which the `series` contains data.
        TODO: Should be a property => end - begin + 1

    """
    def __init__(self, uid, **kwargs):
        super(StationRecord, self).__init__()
        self.uid = uid
        self.ann_anoms = []

    @classmethod
    def invalid(cls, v):
        return v in (MISSING, -MISSING)

    @property
    def series(self):
        """The series of values in celsius."""
        if self._celcius is None:
            c = []
            for v in self._series:
                if self.invalid(v):
                    c.append(XMISSING)
                else:
                    c.append(v * 0.1)
            self._celcius = c
        return self._celcius

    @property
    def series_as_tenths(self):
        """Return the time series in 0.1 celcius, integer units."""
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
        """The shortened form of the record's uinique ID.

        :Return:
            The uid with the country code removed, converted to an integer.

        """
        return int(self.uid[3:])

    def add_year_of_tenths(self, year, data):
        self._add_year_of_data(year, data, MISSING)

    def set_series_from_tenths(self, first_month, series):
        self._set_series(first_month, series, MISSING)

    def set_ann_anoms(self, ann_anoms):
        self.ann_anoms[:] = ann_anoms

    @property
    def ngood_ann_annoms(self):
        """Number of good values in the annual anomolies"""
        # TODO: Rename or remove or something. Also handle magic No.
        if self._ngood_ann_anoms is None:
            bad = 0
            for v in self.ann_anoms:
                bad += v > 9998.99 # TODO: Fix this!
            self._ngood_ann_anoms = len(self.ann_anoms) - bad
        return self._ngood_ann_anoms


    def copy(self):
        r = StationRecord(self.uid)
        r.set_series_from_tenths(self.first_month, self.series_as_tenths)
        return r

    def pretty(self):
        """Format prettily."""
        s = "StationRecord:\n"
        for n, v in sorted(self.__dict__.iteritems()):
            if n.startswith("_"):
                continue
            if n == "meta":
                if v is None:
                    s += "  %-20s: %r\n" % (n, v)
                else:
                    s += "  %-20s: SET\n" % (n, )
            elif n in ("series", "anomolies"):
                s += "  %-20s: [%-4d] = %s...\n" % (n, len(v),
                        str([str(x)[:6] for x in v[:6]])[:-1])
                s += "  %-20s:            ...%s\n" % ("", 
                        str([str(x)[:6] for x in v[-6:]])[1:])
            else:
                s += "  %-20s: %r\n" % (n, v)
        return s[:-1]
           
    # TODO: Yes its a crap name. Currently for debug, but it should be kept
    #       and given a better name.
    def neat(self):
        """Format neatly."""

        s = "StationRecord: %s\n" % (self.station.name)
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
        for n in ("series", "anomolies"):
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

    def pretty(self):
        """Format prettily."""
        s = "SubboxMetaData:\n"
        for n, v in sorted(self.__dict__.iteritems()):
            if not n.startswith("_"):
                s += "  %-20s: %r\n" % (n, v)
        return s[:-1]
           

class SubboxRecord(MonthlyTemperatureRecord):
    """A sub-box record.
    
    This implements the `SubboxProtocol`.

    This is derived from ``tuple``, which provides the iterable requirement
    of the `SubboxProtocol`.

    This can hold, for example, a record of data as stored within the
    ``input/SBBX.HadR2`` file.

    TODO: Describe Ivars.

    Note: Treat `n` as read-only. Use expand to increase it.

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
        #assert self.station_months == self.ngood

    @classmethod
    def invalid(cls, v):
        return abs(v - XMISSING) < 0.1

    @property
    def series(self):
        """The series of values in celsius."""
        return self._series

    def set_series(self, series):
        self._set_series(BASE_YEAR * 12, series, XMISSING)

    @clear_cache
    def pad_with_missing(self, n):
        while self.n < n:
            self._series.append(XMISSING)
            
    @clear_cache
    def set_value(self, idx, value):
        while idx >= len(self.series):
            self._series.append(XMISSING)
        self._series[idx] = value

    # TODO: This seems odd. Do I need this!
    def trim(self):
        self.station_months = len(self._series) - self._series.count(
                XMISSING)

    def __repr__(self):
        return ('<Subbox (%+06.2f,%+06.2f) (%+07.2f,%+07.2f): %d>' %
                (self.lat_S, self.lat_N, self.lon_W, self.lon_E, self.n))
    def neat(self):
        """Format neatly."""

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
        for n in ("series", "anomolies"):
            v = getattr(self, n, None)
            if v is None:
                continue
            s += "   %-10s: [%-4d] = %s...\n" % (n, len(v),
                    str([str(x)[:6] for x in v[:6]])[:-1])
            s += "   %-10s:            ...%s\n" % ("", 
                    str([str(x)[:6] for x in v[-6:]])[1:])

        return s[:-1]

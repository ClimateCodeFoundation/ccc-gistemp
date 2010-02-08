#! /usr/bin/env python
# step0.py
# $URL$
# $Rev$
#
# Nick Barnes, Ravenbrook Limited, 2008-09-17
#
# Copyright (C) 2008-2010 Ravenbrook Limited.  See end of file for license.

"""
Python code for the STEP0 part of the GISTEMP algorithm.

To run:
$ code/step0.py
$

(however, usually run with "python tool/run.py -s 0")

Requires the following files in directory input/:
(tool/fetch.py should know how to get all of these files)

Antarctic/Australasian data:
  antarc1.list
  antarc1.txt
  antarc2.list
  antarc2.txt
  antarc3.list
  antarc3.txt

Hohenpeissenberg station data:
  t_hohenpeissenberg_200306.txt_as_received_July17_2003

USHCN/GHCN station table:
  ushcn2.tbl
  ushcnV2_cmb.tbl

GHCN data:
(from ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/v2/v2.mean.Z)
  v2.mean.Z

USHCN data:
(from ftp://ftp.ncdc.noaa.gov/pub/data/ushcn/v2/monthly/9641C_200907_F52.avg.gz)
  9641C_200907_F52.avg.gz
"""

import itertools
import math
import re

# Should probably be moved elsewhere
def open_or_uncompress(filename):
    """Opens the text file `filename` for reading.  If this fails then
    it attempts to find a compressed version of the file by appending
    '.gz' to the name and opening that (uncompressing it on
    the fly).
    """

    try:
        return open(filename)
    except IOError:
        # When neither  filename nor filename.gz exists we
        # want to pretend that the exception comes from the original
        # call to open, above.  Otherwise the user can be confused by
        # claims that "foo.zip" does not exist when they tried to open
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

# Use the same in-band invalid data marker as step1.  In due course we
# should be able to share code for this sort of thing across all the
# steps (and possibly use a NaN).

def round_value(value):
    """Round a floating-point value to the nearest 0.1. """

    return float(math.floor(value * 10.0 + 0.5)) * 0.1

BAD = round_value(999.9)

def invalid(x):
    """Test for invalid datum ("equal" to the BAD value, for some
    definition of "equal").
    """

    # If you're feeling spooky about the BAD value, re-enable this:
    # if abs(x-BAD) < 0.1:
    #     assert x == BAD

    return x == BAD

def valid(x):
    """Test for valid datum.  See invalid()."""

    # It's important that this obey Iverson's convention: in other words
    # return 0 or 1 (or False or True, which it does).
    return not invalid(x)

### Code to read metadata files.

def read_antarc_station_ids(filename, duplicate):
    """Reads a SCAR station ID file and returns a dictionary mapping
    station name to the 12-digit WMO station ID.  (11 digits from the
    file, plus the final "duplicate flag" digit from our argument).
    """

    dict = {}
    for line in open('input/%s' % filename):
        id11 = line[:11]
        station = line[12:42].strip()
        dict[station] = id11 + duplicate
    return dict

def read_USHCN_stations():
    """Reads the USHCN station list and returns a dictionary
    mapping USHCN station ID to 12-digit WMO station ID.
    """

    stations = {}
    for line in open('input/ushcn2.tbl'):
        (USHCN_id, WMO_id, duplicate) = line.split()
        USHCN_id = int(USHCN_id)
        country_code = int(WMO_id[:3])
        if country_code != 425:
            raise ValueError, "non-425 station found in ushcn.tbl: '%s'" % line
        if duplicate != '0':
            raise ValueError, "station in ushcn.tbl with non-zero duplicate: '%s'" % line
        stations[USHCN_id] = WMO_id+duplicate
    # some USHCNv2 station IDs convert to USHCNv1 station IDs:
    for line in open('input/ushcnV2_cmb.tbl'):
        (v2_station,_,v1_station,_) = line.split()
        stations[int(v2_station)] = stations[int(v1_station)]
    return stations

### Code to read input data files in various formats.

def read_float(s):
    """Converts s to a float, or to BAD if float conversion fails."""

    try:
        return float(s)
    except:
        return BAD

def read_Hohenpeissenberg():
    """reads the Hohenpeissenberg data from
    input/t_hohenpeissenberg_200306.txt_as_received_July17_2003
    which has a header line and then one line per year from 1781.
    We only want data from 1880 to 2002.
    Returns an iterator (id12, year, temperatures).
    """

    for line in open('input/t_hohenpeissenberg_200306.txt_as_received_July17_2003'):
        if line[0] in '12':
            year = int(line[:4])
            if year < 1880 or year > 2002:
                continue
            data = line.split()
            temps = map(read_float, data[1:13])
            temps = temps + [BAD] * (len(temps)-12)
            yield ('617109620002', year, temps)

def read_antarc_line(line):
    """Convert a single line from the Antarctic/Australasian dataset
    files into a year and a 12-tuple of floats (the temperatures in
    Centigrade).  Missing values become BAD.
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

antarc_discard_re = re.compile(r'^$|^Get |^[12A-Z]$')
antarc_temperature_re = re.compile(r'^(.*) .* *temperature')

def read_antarc_data(filename, stations):
    """Reads an Antarctic dataset file - antarc1.txt or antarc3.txt.
    Returns an iterator (id12, year, temperatures).
    """

    for line in open('input/%s' % filename):
        if antarc_discard_re.search(line):
            continue
        station_line = antarc_temperature_re.match(line)
        if station_line:
            station_name = station_line.group(1)
            station_name = station_name.replace('\\','')
            id12 = stations[station_name]
            continue
        line = line.strip()
        if line.find('.') >= 0 and line[0] in '12':
            year, data = read_antarc_line(line)
            yield (id12, year, data)
        
austral_discard_re = re.compile(r'^$|:')
austral_header_re = re.compile(r'^\s*(.+?)  .*(E$|E )')

def read_austral_data(filename, stations):
    """Reads the Australasian dataset file antarc2.txt.
    Returns an iterator (id12, year, temperatures).
    """

    for line in open('input/%s' % filename):
        if austral_discard_re.search(line):
            continue
        station_line = austral_header_re.match(line)
        if station_line:
            station_name = station_line.group(1).strip()
            id12 = stations[station_name]
            continue
        line = line.strip()
        if line.find('.') >= 0 and line[0] in '12':
            year, data = read_antarc_line(line)
            yield (id12, year, data)
        
def read_antarctica():
    """Reads the SCAR datasets (Antarctic and Australasian).
    Returns an iterator (id12, year, temperatures).
    """

    stations1 = read_antarc_station_ids('antarc1.list', '8')
    stations2 = read_antarc_station_ids('antarc2.list', '7')
    stations3 = read_antarc_station_ids('antarc3.list', '9')
    return itertools.chain(read_antarc_data('antarc1.txt', stations1),
                           read_antarc_data('antarc3.txt', stations3),
                           read_austral_data('antarc2.txt', stations2))

def read_GHCN(filename):
    """Reads a GHCN data file (usually v2.mean).
    Returns an iterator (id12, year, temperatures).
    """

    GHCN_last_year = 0
    for line in open_or_uncompress('input/%s' % filename):
        id12 = line[:12]
        year = int(line[12:16])
        GHCN_last_year = max(year, GHCN_last_year)
        temps = []
        for i in range(0,12):
            t = float(line[i*5+16:i*5+21])
            if t < -9998:
                temps.append(BAD)
            else:
                temps.append(t / 10.0)
        yield (id12, year, temps)
    t = open('work/GHCN.last_year','w')
    t.write('%d\n' % GHCN_last_year)
    t.close()

def read_USHCN():
    """Returns USHCN dataset as a dictionary mapping 12-digit WMO
    station ID to a dictionary mapping year number to temperature data
    list.
    """

    hash = {}
    last_year = 0
    stations = read_USHCN_stations()
    # initialize hash from WMO station ID to results
    for (USHCN_id, id12) in stations.items():
        hash[id12] = {}
    for line in open_or_uncompress('input/9641C_200907_F52.avg'):
        if line[6] != '3': # 3 indicates mean temperatures
            continue
        USHCN_station = int(line[0:6])
        id12 = stations[USHCN_station]
        year = int(line[7:11])
        if year < 1880: # discard data before 1880
            continue
        last_year = max(year, last_year)
        temps = []
        valid = False
        for m in range(0,12):
            temp_fahrenheit = int(line[m*7+11:m*7+17])
            flag = line[m*7+17]
            if ((flag in 'EQ') or              # interpolated data
                (temp_fahrenheit == -9999)) :  # absent data
                temp = BAD
            else:
                # degrees centigrade
                temp = (temp_fahrenheit - 320) * 5/90.0
                valid = True
            temps.append(temp)
        if valid: # some valid data found
            hash[id12][year] = temps
    return (hash, last_year)

### Code to combine and manipulate input data sets

def include_US(GHCN):
    """Adds the USHCN data to the GHCN dataset (passed in as an
    iterator argument).  Returns an iterator for the combined dataset.
    Where the GHCN dataset contains data for stations in the USHCN
    set, return the USHCN data, adjusted for differences in
    monthly temperature means.
    """

    GHCN_hash = {}
    (USHCN, last_year) = read_USHCN()
    while True:
        try:
            GHCN_line = GHCN.next()
            (id12, year, data) = GHCN_line
            if USHCN.has_key(id12):
                if not USHCN[id12].has_key(year):
                    # discard this GHCN result line
                    continue
                GHCN_hash.setdefault(id12,{})
                GHCN_hash[id12][year] = data
            else:
                yield GHCN_line
        except StopIteration: # end of GHCN iterator
            # now iterate over USHCN stations.
            for (id12, udata) in USHCN.items():
                gdata = GHCN_hash.get(id12, {})
                # Calculate a set of monthly offsets USHCN-GHCN for
                # this station.
                # For compatibility with GISTEMP, these offsets are
                # averages over up to ten years, working backwards
                # from the present day to 1979, and are rounded to the
                # nearest tenth of a degree.
                diffs = []
                for month in range(0,12):
                    diff = 0
                    count = 0
                    for year in range(last_year, 1979, -1):
                        if (udata.has_key(year) and valid(udata[year][month]) and
                            gdata.has_key(year) and valid(gdata[year][month])):
                            diff += udata[year][month]-gdata[year][month]
                            count += 1
                            if count == 10:
                                break
                    if count > 0:
                        # compatibility with GISTEMP requires that we
                        # round this to the nearest tenth of a degree.
                        average = diff/count
                        diffs.append(average)
                    else:
                        diffs.append(0)
                # Now apply the USHCN-GHCN offsets to every year for
                # this station in which there is GHCN data.  Note that
                # years in USHCN but not in GHCN are just copied, not
                # adjusted.  This may be a bug.
                for year in sorted(udata.keys()):
                    temps = udata[year]
                    if gdata.has_key(year):
                        # station-years in GHCN are adjusted
                        def adjust(t,d):
                            if valid(t):
                                return t-d
                            else:
                                return BAD
                        adjusted = map(adjust, temps, diffs)
                        yield(id12, year, adjusted)
                    else:
                        yield(id12, year, temps)
            return

def dump_old(set, year):
    """Returns the input set without any items in years prior to the
    specified year.
    """
    
    return (item for item in set if item[1] >= year)

def remove_Hohenpeissenberg(set):
    """Remove data from an iterator if it from Hohenpeisenberg prior
    to 2003.
    """

    return (t for t in set if (t[0][:11] != '61710962000' or
                               t[1] > 2002))

def from_years(years):
    """*years* is a list of year records (lists of temperatures) that
    comprise a station's entire record.  The data are converted to a
    linear array (could be a list/tuple/array/sequence, I'm not
    saying), *series*, where series[m] gives the temperature (a
    floating point value in degrees C) for month *m*, counting from 0
    for the January of the first year with data.

    (*series*,*begin*) is returned, where *begin* is
    the first year for which there is data for the station.

    Code almost identical to this is also in step1.py at present, and
    should be shared.
    """

    begin = None
    # Previous year.
    prev = None
    series = []
    for (year, data) in years:
        if begin is None:
            begin = year
        # Some input datasets can have duplicate years.  We just drop
        # the duplicate data here.
        if prev and prev >= year:
            continue
        # The sequence of years for a station record is not
        # necessarily contiguous.  For example "1486284000001988" is
        # immediately followed by "1486284000001990", missing out 1989.
        # Extend with blank data.
        while prev and prev < year-1:
            series.extend([BAD]*12)
            prev += 1
        prev = year
        series.extend(data)
    return (series, begin)

def as_station_record(id, lines):
    """Takes a set of *lines*, which are triples (id, year, data), for
    a single station, and returns a pair *(dict, series)*, in which
    *series* is a combined list of monthly values and *dict* is a
    metadata dictionary.
    """

    (series, begin) = from_years((year, data) for (_, year, data) in lines)
    dict = {}
    dict['id'] = id
    dict['begin'] = begin
    return (dict, series)


### Main step 0 function.

def step0():
    """Replaces GISTEMP STEP0.
    Reads the GHCN, USHCN, SCAR, and Hohenpeissenberg data
    and returns a single combined dataset as an iterator
    of *(metadata, series)* pairs.
    """

    combined = itertools.chain(read_GHCN("v2.mean"),
                               read_antarctica())
    since_1880 = dump_old(combined, 1880)
    including_US = include_US(since_1880)
    with_Hohenpeissenberg = itertools.chain(read_Hohenpeissenberg(),
                                            remove_Hohenpeissenberg(including_US))
    sorted_lines = sorted(with_Hohenpeissenberg)
    for (id, lines) in itertools.groupby(sorted_lines, lambda line: line[0]):
        yield as_station_record(id, lines)

# Code to output result set; these functions will move out of step0.py
# in due course.

def round_to_nearest(f):
    """Returns the int which is nearest to the argument float.  Draws
    are rounded away from zero.
    """

    if f >= 0:
        return int(math.floor(f + 0.5))
    else:
        return int(math.ceil(f - 0.5))

def str_temp(temp):
    """Takes a float, and converts to 5 characters of ASCII: decimal
    representation or "-9999" if None.
    """

    if invalid(temp):
        return '-9999'
    else:
        return "%5d" % round_to_nearest(temp * 10.0)

def str_line(id, year, temps):
    """Turns a single year from the dataset into a line in the v2.mean format."""

    return ('%s%04d%s\n' % (id, year, ''.join(map(str_temp, temps))))

def write_data(set, filename):
    """Writes a dataset to a file in v2.mean format."""

    f = open(filename, 'w')
    for (dict, series) in set:
        id = dict['id']
        year = dict['begin']
        for months in (series[i:i+12] for i in range(0,len(series),12)):
            f.write(str_line(id, year, months))
            year = year + 1
    f.close()

def main():
    write_data(step0(), "work/v2.mean_comb")

if __name__ == '__main__':
    main()

# This file is copyright (C) 2008 Ravenbrook Limited.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1.  Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#
# 2.  Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in
#     the documentation and/or other materials provided with the
#     distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDERS AND CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
# TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
# USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.

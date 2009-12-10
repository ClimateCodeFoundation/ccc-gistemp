#! /usr/bin/env python
#
# step0.py
# $Id: //info.ravenbrook.com/project/ccc/master/code/step0.py#10 $
# Copyright (C) 2008 Ravenbrook Limited.  See end of file for license.
# 
# Python code for the STEP0 part of the GISTEMP algorithm.
# Nick Barnes, Ravenbrook Limited, 2008-09-17
#
# To run:
# $ code/step0.py
# $ sort < work/v2.mean_comb.unsorted > work/v2.mean_comb
# $
#
# Requires the following files in directory input/:
#
# From GISTEMP/STEP0/input_files, Antarctic/Australasian data:
#   antarc1.list
#   antarc1.txt
#   antarc2.list
#   antarc2.txt
#   antarc3.list
#   antarc3.txt
#
# From GISTEMP/STEP0/input_files, Hohenpeissenberg station data:
#   t_hohenpeissenberg_200306.txt_as_received_July17_2003
#
# From GISTEMP/STEP0/input_files, USHCN/GHCN station table:
#   ushcn.tbl
#
# From import/2008-07-14/ghcn, GHCN data:
# (originally from ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/v2/v2.mean.Z)
#   v2.mean
#
# From import/2007-10-11/ushcn, USHCN data:
# (originally from ftp://ftp.ncdc.noaa.gov/pub/data/ushcn/hcn_doe_mean_data.Z)
#   hcn_doe_mean_data

import itertools
import math
import re

# Should probably be moved elsewhere
def open_or_uncompress(filename):
    """Opens the text file `filename` for reading.  If this fails then
    it attempts to find a compressed version of the file by appending
    '.Z' to the name and opening that (uncompressing it on the fly).
    """

    try:
        return open(filename)
    except IOError:
        # When neither filename nor filename.Z exists we want to pretend
        # that the exception comes from the original call to open,
        # above.  Otherwise the user can be confused by claims that
        # "foo.Z" does not exist when they tried to open "foo".  In
        # order to maintain this pretence, we have to get the exception
        # info and save it. See
        # http://blog.ianbicking.org/2007/09/12/re-raising-exceptions/
        import sys
        exception = sys.exc_info()
        try:
            import uncompress
            return uncompress.Zfile(filename + '.Z')
        except IOError:
            pass
        raise exception[0], exception[1], exception[2]

def read_antarc_station_ids(filename):
    """Reads a SCAR station ID file and returns a dictionary
    mapping station name to the WMO triple
    (country_code, WMO_station, modifier).
    """

    dict = {}
    for line in open('input/%s' % filename):
        country_code = int(line[:3])
        WMO_station = int(line[3:8])
        modifier = int(line[8:11])
        station = line[12:42].strip()
        dict[station] = (country_code, WMO_station, modifier)
    return dict

def read_float(s):
    """Converts s to a float, or to None if float conversion fails."""

    try:
        return float(s)
    except:
        return None

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

def read_tenths(s):
    """Returns the integer nearest to the argument string times 10.
    If float conversion fails, returns None.
    """

    try:
        f = float(s)
    except:
        return None
    return round_to_nearest(f * 10)

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

antarc_discard_re = re.compile(r'^$|^Get |^[12A-Z]$')
antarc_temperature_re = re.compile(r'^(.*) .* *temperature')

def read_antarc_data(filename, duplicate, stations):
    """An iterator to read an Antarctic dataset file -
    antarc1.txt or antarc3.txt.
    """

    for line in open('input/%s' % filename):
        if antarc_discard_re.search(line):
            continue
        station_line = antarc_temperature_re.match(line)
        if station_line:
            station_name = station_line.group(1)
            station_name = station_name.replace('\\','')
            (country_code, WMO_station, modifier) = stations[station_name]
            continue
        line = line.strip()
        if line.find('.') >= 0 and line[0] in '12':
            year, data = read_antarc_line(line)
            yield (country_code, WMO_station, modifier, duplicate, year, data)
        
austral_discard_re = re.compile(r'^$|:')
austral_header_re = re.compile(r'^\s*(.+?)  .*(E$|E )')

def read_austral_data(filename, duplicate, stations):
    """An iterator to read the Australasian dataset file antarc2.txt."""

    for line in open('input/%s' % filename):
        if austral_discard_re.search(line):
            continue
        station_line = austral_header_re.match(line)
        if station_line:
            station_name = station_line.group(1).strip()
            (country_code, WMO_station, modifier) = stations[station_name]
            continue
        line = line.strip()
        if line.find('.') >= 0 and line[0] in '12':
            year, data = read_antarc_line(line)
            yield (country_code, WMO_station, modifier, duplicate, year, data)
        
def read_antarctica():
    """An iterator for the SCAR dataset (Antarctic and Australasian)."""

    stations1 = read_antarc_station_ids('antarc1.list')
    stations2 = read_antarc_station_ids('antarc2.list')
    stations3 = read_antarc_station_ids('antarc3.list')
    return itertools.chain(read_antarc_data('antarc1.txt', 8, stations1),
                           read_antarc_data('antarc3.txt', 9, stations3),
                           read_austral_data('antarc2.txt', 7, stations2))

def read_GHCN(filename):
    """An iterator to read a GHCN data file (usually v2.mean)."""

    GHCN_last_year = 0
    for line in open_or_uncompress('input/%s' % filename):
        country_code = int(line[:3])
        WMO_station = int(line[3:8])
        modifier = int(line[8:11])
        duplicate = int(line[11:12])
        year = int(line[12:16])
        GHCN_last_year = max(year, GHCN_last_year)
        temps = []
        for i in range(0,12):
            t = int(line[i*5+16:i*5+21])
            if t == -9999:
                temps.append(None)
            else:
                temps.append(t)
        yield((country_code, WMO_station, modifier, duplicate, year, temps))
    t = open('work/GHCN.last_year','w')
    t.write('%d\n' % GHCN_last_year)

def dump_old(set, year):
    """Returns the input set without any items in years prior to the
    specified year.
    """
    
    return (item for item in set if item[4] >= year)

def read_USHCN_stations():
    """Reads the USHCN station list and returns a dictionary
    mapping USHCN station ID to (WMO station ID, modifier).
    """

    stations = {}
    for line in open('input/ushcn.tbl'):
        # USHCN2v2.f: read(1,'(5x,i6,1x,a11,1x,i1)') idus,idg,iz
        USHCN_id = int(line[5:11])
        country_code = int(line[12:15])
        WMO_station = int(line[15:20])
        modifier = int(line[20:23])
        duplicate = int(line[25])
        if country_code != 425:
            raise ValueError, "non-425 station found in ushcn.tbl: '%s'" % line
        if duplicate != 0:
            raise ValueError, "station in ushcn.tbl with non-zero duplicate: '%s'" % line
        stations[USHCN_id] = (WMO_station, modifier)
    return stations

def read_USHCN():
    """Returns USHCN dataset as a dictionary mapping WMO station ID to
    a dictionary mapping year number to temperature data list.
    """

    hash = {}
    last_year = 0
    stations = read_USHCN_stations()
    # initialize hash from WMO station ID to results
    for (USHCN_id, (WMO_station, modifier)) in stations.items():
        hash[(WMO_station, modifier)] = {}
    for line in open_or_uncompress('input/hcn_doe_mean_data'):
        # " 3A" indicates mean Filnet temperatures.
        if line[11:14] != ' 3A':
            continue
        USHCN_station = int(line[0:6])
        (WMO_station, modifier) = stations[USHCN_station]
        year = int(line[7:11])
        if year < 1880: # discard data before 1880
            continue
        last_year = max(year, last_year)
        temps = []
        valid = False
        for m in range(0,12):
            temp_fahrenheit = read_float(line[m*10+14:m*10+20])
            flags = line[m*10+20:m*10+24]
            if ((flags[3] == 'M') or           # interpolated data
                (temp_fahrenheit < -99)) :     # absent data
                temp = None
            else:
                # tenths of degree centigrade
                temp = round_to_nearest((temp_fahrenheit - 32) * 50/9.0)
                valid = True
            temps.append(temp)
        if valid: # some valid data found
            hash[(WMO_station, modifier)][year] = temps
    return (hash, last_year)

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
            (country_code, WMO_station, modifier, duplicate, year, data) = GHCN_line
            key = (WMO_station, modifier)
            if country_code == 425 and duplicate == 0 and USHCN.has_key(key):
                if not USHCN[key].has_key(year):
                    # discard this GHCN result line
                    continue
                GHCN_hash.setdefault(key,{})
                GHCN_hash[key][year] = data
            else:
                yield GHCN_line
        except: # end of GHCN iterator
            # now iterate over USHCN stations.
            for (key, udata) in USHCN.items():
                (WMO_station, modifier) = key
                gdata = GHCN_hash.get(key, {})
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
                        if (udata.has_key(year) and udata[year][month] is not None and
                            gdata.has_key(year) and gdata[year][month] is not None):
                            diff += udata[year][month]-gdata[year][month]
                            count += 1
                            if count == 10:
                                break
                    if count > 0:
                        # compatibility with GISTEMP requires that we
                        # round this to the nearest tenth of a degree.
                        average = round_to_nearest(float(diff)/count)
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
                            if t is None:
                                return None
                            else:
                                return t-d
                        adjusted = map(adjust, temps, diffs)
                        yield(425, WMO_station, modifier, 0, year, adjusted)
                    else:
                        yield(425, WMO_station, modifier, 0, year, temps)
            return

def read_Hohenpeissenberg():
    """reads the Hohenpeissenberg data from
    input/t_hohenpeissenberg_200306.txt_as_received_July17_2003
    which has a header line and then one line per year from 1781.
    We only want data from 1880 to 2002.
    """

    for line in open('input/t_hohenpeissenberg_200306.txt_as_received_July17_2003'):
        if line[0] in '12':
            year = int(line[:4])
            if year < 1880 or year > 2002:
                continue
            data = line.split()
            temps = map(read_tenths, data[1:13])
            temps = temps + [None] * (len(temps)-12)
            yield (617, 10962, 0, 2, year, temps)

def remove_Hohenpeissenberg(set):
    """Remove data from an iterator if it from Hohenpeisenberg prior
    to 2003.
    """

    return (t for t in set if (t[0] != 617 or
                               t[1] != 10962 or
                               t[2] != 0 or
                               t[4] > 2002))

def str_temp(temp):
    """Takes either None or an integer, and converts to 5 characters
    of ASCII: decimal representation or "-9999" if None.
    """

    if temp is None:
        return '-9999'
    else:
        return "%5d" % temp

def str_line(line):
    """Turns a dataset item into a line in the v2.mean format."""

    (country_code, WMO_station, modifier, duplicate, year, temps) = line
    return ('%03d%05d%03d%1d%04d%s\n' %
            (country_code, WMO_station, modifier,
             duplicate, year, ''.join(map(str_temp, temps))))

def write_data(set, filename):
    """Writes a dataset to a file in v2.mean format."""

    f = open(filename, 'w')
    for line in set:
        f.write(str_line(line))
    f.close()

def step0(input_file):
    """Replaces GISTEMP STEP0.
    Reads the GHCN, USHCN, SCAR, and Hohenpeissenberg data
    and returns a single combined dataset as an iterator.
    """

    combined = itertools.chain(read_GHCN(input_file),
                               read_antarctica())
    since_1880 = dump_old(combined, 1880)
    including_US = include_US(since_1880)
    return itertools.chain(read_Hohenpeissenberg(),
                           remove_Hohenpeissenberg(including_US))

def main():
    write_data(step0("v2.mean"),
               "work/v2.mean_comb.unsorted")

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

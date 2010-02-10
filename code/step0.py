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

import giss_data

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


def read_GHCN(ghcn_source, antarc_source):
    """An iterator to read a GHCN data file (usually v2.mean)."""
    for record in ghcn_source:
        yield record
    for record in antarc_source:
        yield record


def calc_monthy_USHCN_offsets(u_record, g_record):
    u_years = u_record.get_set_of_years_as_tenths(1980, u_record.last_year)
    g_years = g_record.get_set_of_years_as_tenths(1980, u_record.last_year)
    reversed_year_pairs = list(reversed(zip(u_years, g_years)))

    diffs = [0] * 12
    for month in range(12):
        sum = 0.0
        count = 0
        for u_year, g_year in reversed_year_pairs:
            u_temp, g_temp = u_year[month], g_year[month]
            if giss_data.MISSING not in (u_temp, g_temp):
                sum += u_temp - g_temp
                count += 1
                if count == 10:
                    break

        if count > 0:
            # compatibility with GISTEMP requires that we round this
            # to the nearest tenth of a degree.
            diffs[month] = round_to_nearest(sum / count)

    return diffs


def adjust_USHCN(ushcn_records, ghcn_records, us_only):
    """Where the GHCN dataset contains data for stations in the USHCN
    set, adjust the USHCN data, for differences in
    monthly temperature means.
    """
    print "Adjust USHCN records"

    def adj(t, d):
        if t != giss_data.MISSING:
            return t - d
        return t

    # For each USHCN record look for a corresponding GHCN record. Matching
    # GHCN records have the same station_uid and a discriminator of zero.
    # Hence we can generate a lookup key by: ``u_record.station_uid + "0"``.
    to_remove = []
    for u_key, u_record in ushcn_records.iteritems():
        key = u_record.station_uid + "0"
        g_record = ghcn_records.get(key, None)
        if g_record is None:
            print "NO MODS", u_record.uid
            continue
        diffs = calc_monthy_USHCN_offsets(u_record, g_record)

        # Now apply the USHCN-GHCN offsets to every year for
        # this station in which there is GHCN data.  Note that
        # years in USHCN but not in GHCN are just copied, not
        # adjusted.  This may be a bug.
        tenths = list(u_record.series_as_tenths)
        new_data = []
        for year in range(u_record.first_year, u_record.last_year + 1):
            temps = u_record.get_a_year_as_tenths(year)
            if g_record.has_data_for_year(year):
                new_data.extend([adj(t, d) for t, d in zip(temps, diffs)])
            else:
                new_data.extend(temps)
        g_record.set_series_from_tenths(u_record.first_year * 12 + 1, new_data)
        to_remove.append(u_key)
        us_only.pop(key, None)

    for k in to_remove:
        del ushcn_records[k]
    for k in us_only:
        ghcn_records.pop(k, None)


def remove_Hohenpeissenberg_from_GHCN(ghcn_records, hohenpeis_record):
    print "Remove Hohenpeissenberg data from GHCN records"

    for g_record in ghcn_records.itervalues():
        if g_record.station_uid == hohenpeis_record.station_uid:
            # Extract the data for the years 2003 to present.
            new_data = []
            for year in range(2003, g_record.last_year + 1):
                new_data.extend(g_record.get_a_year_as_tenths(year))

            if g_record.uid == hohenpeis_record.uid:
                # If the record has the same UID as Hohenpeissenberg then
                # replace the data with Hohenpeissenberg plus the recent
                # years.
                last_year = g_record.last_year
                g_record.set_series_from_tenths(hohenpeis_record.first_month,
                        hohenpeis_record.series_as_tenths)
                for i, year in enumerate(range(2003, last_year + 1)):
                    g_record.add_year_of_tenths(year, new_data[i * 12:(i + 1) * 12])

            else:
                # For other records, just replace the data with the later
                # years.
                g_record.set_series_from_tenths(2002 * 12 + 1, new_data)


# TODO: Should the antarc_source really appear at this stage?
def asdict(ghcn_source, antarc_source):
    print "Load GHCN records"
    return dict((record.uid, record) for record in
      itertools.chain(ghcn_source, antarc_source))


def read_USHCN(ushcn_source):
    """Read all USHCN station records from the input.

    :Return:
        A tuple of ``(records, us_only)``. The ``records`` is a dictionary of
        all the records, keyed by the record UID. The ``us_only`` is a
        dictionary, where each key is the UID that would match a a
        corresponding station record from the GHCN data source.

    """
    print "Load USHCN records"
    d = dict((record.uid, record) for record in ushcn_source)
    return d, ushcn_source.us_only


def step0(inputs):
    """An iterator for step 1.  Produces a stream of
    `giss_data.StationRecord` instances.

    """
    ushcn_records, us_only = read_USHCN(inputs.ushcn_source)
    ghcn_records = asdict(inputs.ghcn_source, 
            inputs.antarc_source)
    adjust_USHCN(ushcn_records, ghcn_records, us_only)
    hohenpeis_record = list(inputs.hohenpeis_source)[0]
    remove_Hohenpeissenberg_from_GHCN(ghcn_records, hohenpeis_record)

    records = {}
    records.update(ghcn_records)
    records.update(ushcn_records)

    # TODO: I think the whole point of sorting here is so that step1
    # receives the records in a good order to perfrom its
    # ``comb_records(self.record_source)`` processing. We can probably
    # dispense with this by making use of the fact that we have a
    # dictionary here, so we can combine at this stage. [Paul O]
    for uid, record in sorted(records.iteritems()):
        if record.n:
            yield record


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

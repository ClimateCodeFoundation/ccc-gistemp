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
$ python tool/run.py -s 0

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
(from ftp://ftp.ncdc.noaa.gov/pub/data/ushcn/v2/monthly/9641C_*_F52.avg.gz)
  ushcnv2.gz
"""

# http://docs.python.org/release/2.4.4/lib/module-itertools.html
import itertools

# Clear Climate Code
from giss_data import valid
import parameters

def calc_monthly_USHCN_offsets(u_record, g_record):
    u_years = u_record.get_set_of_years(parameters.USHCN_offset_start_year,
                                        u_record.last_year)
    g_years = g_record.get_set_of_years(parameters.USHCN_offset_start_year,
                                        u_record.last_year)
    reversed_year_pairs = list(reversed(zip(u_years, g_years)))

    diffs = [0.0] * 12
    for month in range(12):
        sum = 0.0
        count = 0
        for u_year, g_year in reversed_year_pairs:
            u_temp, g_temp = u_year[month], g_year[month]
            if valid(u_temp) and valid(g_temp):
                sum += u_temp - g_temp
                count += 1
                if count == parameters.USHCN_offset_max_months:
                    break
        if count > 0:
            diffs[month] = sum / count
    return diffs


def adjust_USHCN(ushcn_records, ghcn_records, us_stations):
    """Where the GHCN dataset, *ghcn_records*, contains data for
    stations in the USHCN dataset, *ushcn_records*, replace the
    GHCN data with USHCN data, adjusted for differences in monthly
    temperature means.  The GHCN data is updated; the corresponding
    USHCN station is removed from *ushcn_records*.  *us_stations*
    should be a collection of 12-digit GHCN record identifiers,
    normally these are the records that are replaced by USHCN data,
    but if they are not replaced by USHCN data, they are removed
    (from *ghcn_records*).
    """

    print "Adjust USHCN records"

    us_only = {}
    for id12 in us_stations:
        us_only[id12] = None

    def adj(t, d):
        if valid(t):
            return t - d
        return t

    # For each USHCN record look for a corresponding GHCN record.
    to_remove = []
    for key, u_record in ushcn_records.iteritems():
        g_record = ghcn_records.get(key, None)
        if g_record is None:
            print "NO MODS", u_record.uid
            continue
        diffs = calc_monthly_USHCN_offsets(u_record, g_record)

        # Now apply the USHCN-GHCN offsets to every year for
        # this station in which there is GHCN data.  Note that
        # years in USHCN but not in GHCN are just copied, not
        # adjusted.  This may be a bug.
        new_data = []
        for year in range(u_record.first_year, u_record.last_year + 1):
            temps = u_record.get_a_year(year)
            if g_record.has_data_for_year(year):
                new_data.extend([adj(t, d) for t, d in zip(temps, diffs)])
            else:
                new_data.extend(temps)
        g_record.set_series(u_record.first_year * 12 + 1, new_data)
        to_remove.append(key)
        us_only.pop(key, None)

    for k in to_remove:
        del ushcn_records[k]
    for k in us_only:
        ghcn_records.pop(k, None)

def correct_Hohenpeissenberg(ghcn_records, hohenpeissenberg):
    """Replace Hohenpeissenberg data from 1880 to 2002 in the GHCN
    dataset with the priv. comm. data."""
    print "Correct the GHCN Hohenpeissenberg record."

    cut = hohenpeissenberg.last_year + 1

    for record in ghcn_records.itervalues():
        if record.station_uid == hohenpeissenberg.station_uid:
            # Extract the data for the years more recent than the priv.
            # comm. data.
            new_data = []
            for year in range(cut, record.last_year + 1):
                new_data.extend(record.get_a_year(year))

            if record.uid == hohenpeissenberg.uid:
                # If the record has the same UID as Hohenpeissenberg then
                # replace the data with Hohenpeissenberg plus the recent
                # years.
                last_year = record.last_year
                record.set_series(hohenpeissenberg.first_month,
                                  hohenpeissenberg.series)
                for i, year in enumerate(range(cut, last_year + 1)):
                    record.add_year(year, new_data[i * 12:(i + 1) * 12])

            else:
                # For other records, just replace the data with the later
                # years.
                record.set_series(cut * 12 + 1, new_data)


def asdict(records):
    """Convert simple series of stations records into a dictionary
    (mapping from uid to record)."""

    return dict((record.uid, record) for record in records)


def step0(inputs):
    """An iterator for step 0.  Produces a stream of
    `giss_data.Series` instances.

    """
    ushcn_records = asdict(inputs.ushcn_source)
    print "Load GHCN records"
    ghcn_records = asdict(inputs.ghcn_source)
    scar_records = asdict(inputs.antarc_source)
    if 'hohenpeissenberg' in parameters.data_sources:
        correct_Hohenpeissenberg(ghcn_records, inputs.hohenpeissenberg)
    if ('ushcn' in parameters.data_sources and
        'ghcn' in parameters.data_sources):
        adjust_USHCN(ushcn_records, ghcn_records,
          inputs.ushcn_stations.itervalues())

    records = {}
    if 'ghcn' in parameters.data_sources:
        records.update(ghcn_records)
    if 'ushcn' in parameters.data_sources:
        records.update(ushcn_records)
    if 'scar' in parameters.data_sources:
        records.update(scar_records)

    # We sort here - as does GISTEMP - so that all the records for a
    # given 11-digit station ID are grouped together, ready for
    # combining in the next step.  I believe our inputs already have
    # this property (and we can assert it there as we go along), so we
    # can probably ensure it here by doing the above processing with
    # suitable care.
    for uid, record in sorted(records.iteritems()):
        if record:
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

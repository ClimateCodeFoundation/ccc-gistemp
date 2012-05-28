#! /usr/bin/env python
# $URL$
# $Rev$
#
# extension/step0.py
#
# Nick Barnes, Climate Code Foundation, 2012-05-27
#
# Copyright (C) Ravenbrook Limited, 2008-2010.
# Copyright (C) Climate Code Foundation, 2010-2012.
# BSD license, see license.txt

"""
Python code for extensions to the Step 0 part of the GISTEMP
algorithm.  This includes functionality to combine the USHCN and GHCN
datasets, which were used by GISTEMP before they were both
superceded in 2011-12 by GHCN-M version 3.
"""
import os

import parameters
from code import step0
from code.giss_data import valid

log = open(os.path.join('log', 'step0.log'), 'w')

def calc_monthly_USHCN_offsets(u_record, g_record):
    """Given a USHCN record `u_record` and a GHCN record `g_record`,
    using any overlapping years, computes a set of 12 monthly offsets
    which can be added to the GHCN record to most closely approximate
    the USHCN record."""
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


def adjust_USHCN(ushcn_records, ghcn_records):
    """Each USHCN record is adjusted for differences in monthly means
    between it and the corresponding record in GHCN.  The corresponding
    GHCN record is removed.  USHCN records that have no
    corresponding GHCN records are retained unadjusted (ordinarily there
    are no such records, but may be if subsets are created by hand).
    """

    def adj(t, d):
        if valid(t):
            return t - d
        return t

    # Count of USHCN records that have no GHCN counterpart.
    count_ushcn_only = 0
    # For each USHCN record look for a corresponding GHCN record.
    for key, u_record in ushcn_records.iteritems():
        g_record = ghcn_records.get(key, None)
        if g_record is None:
            count_ushcn_only += 1
            log.write("""%s action "ushcn only"\n""" % key)
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
        u_record.set_series(u_record.first_year * 12 + 1, new_data)
        del ghcn_records[key]

    if count_ushcn_only:
        print count_ushcn_only, "USHCN records had no GHCN counterpart."

def discard_contig_us(records):
    """Discard stations in the contiguous US.  The stations are
    discarded on the basis of their 11-digit GHCN ID: stations from
    425710000000 to 425900000000 are excluded, not including the latter
    end point.  Note: US stations below this range are in Alaska, US
    stations after this range are Pacific Ocean islands.
    """

    for uid in records.keys():
        assert len(uid) == 12
        if '425710000000' <= uid < '425900000000':
            # Note: records.keys() produces a fresh list and that's
            # important for making this del safe.
            del records[uid]

def pre_step0(input):
    """Apply whatever extensions we have for GISTEMP step 0, that run
    before the main step 0.  At present this is used for two purposes.
    Firstly: incorporating the GHCN version 2 dataset, when it is
    specified in the parameters, with the same semantics as GHCN-M
    version 3 Secondly: managing the interaction between USHCN and
    GHCN, if both datasets are present, as was normal for GISTEMP
    before 2011-12.
    """

    real_sources = input.sources
    if 'ghcn.v2' in input.sources:
        print "Extension: use GHCN version 2."
        input.sources = real_sources.copy()
        input.sources[input.sources.index('ghcn.v2')] = 'ghcn'
        ghcn_key = 'ghcn.v2'
    else:
        ghcn_key = 'ghcn'

    real_open = input.open
    if 'ghcn' in input.sources and 'ushcn' in input.sources:
        print "Extension: merge USHCN and GHCN records."
        # Read GHCN and USHCN datasets into big temporary lists
        ghcn_records = list(input.open(ghcn_key))
        ushcn_records = list(input.open('ushcn'))
        # Turn the lists into dictionaries for processing
        ghcn_data = dict((record.uid, record) for record in ghcn_records)
        ushcn_data = dict((record.uid, record) for record in ushcn_records)
        # Remember the ordering of the records ...
        ghcn_ids = [record.uid for record in ghcn_records]
        ushcn_ids = [record.uid for record in ushcn_records]
        # ... but forget the big lists themselves, as they
        # can be very large.
        ghcn_records = None
        ushcn_records = None

        # function to turn (list of ids, dict) back into records
        # generator
        def data(ids, dict):
            for id in ids:
                if id in dict:
                    yield dict[id]

        # ... adjust the USHCN data;
        # (this can affect both dicts)
        adjust_USHCN(ushcn_data, ghcn_data)

        # ... optionally exclude contiguous US stations.
        if parameters.retain_contiguous_US:
            print "Extension: retain US data in GHCN."
        else:
            discard_contig_us(ghcn_data)

        # build a fake open method
        # that gives back these modified datasets
        # (and nukes the temporaries).
        def fake_open(source):
            if source == 'ghcn':
                return data(ghcn_ids, ghcn_data)
                ghcn_ids = None
                ghcn_data = None
            elif source == 'ushcn':
                return data(ushcn_ids, ushcn_data)
                ushcn_ids = None
                ushcn_data = None
            else:
                return real_open(source)

    else:
        def fake_open(source):
            if source == 'ghcn':
                return real_open(ghcn_key)
            else:
                return real_open(source)

    input.open = fake_open

    return input

def post_step0(records):
    """Apply whatever extensions we have for GISTEMP step 0, that run
    after the main step 0.  None at present."""
    return records

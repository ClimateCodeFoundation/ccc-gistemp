#! /usr/bin/env python
# $URL$
# $Rev$
#
# step1.py
#
# Nick Barnes, Ravenbrook Limited, 2008-08-06
# David Jones, Climate Code Foundation, 2010-09-22

"""
Python code reproducing the Step 1 part of the GISTEMP algorithm.  In
this step: duplicate records from the same station (a GHCN v2 feature)
are combined into a single record, if possible; certain records (or
parts) are adjusted or dropped, under the control of configuration
files.

Requires the following files in the input/ directory,
from GISTEMP STEP1/input_files/:

mcdw.tbl
ushcn2.tbl
sumofday.tbl
v2.inv

Do not be tempted to replace v2.inv with the apparently similar
v2.temperature.inv file available from NOAA,
ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/v2/v2.temperature.inv .  The
NOAA file has been treated for GISTEMP's use by, for example, adding
records corresponding to Antarctic stations that are not used in GHCN
but are used in the GISTEMP analysis.  Step 1 (this step) expects to
find a record in v2.inv for every station it has a time series for.

Requires the following files in the config/ directory,
from GISTEMP STEP1/input_files/:

combine_pieces_helena.in
Ts.strange.RSU.list.IN
Ts.discont.RS.alter.IN

Also requires the existence of writeable work/ and log/ directories.
"""

import math
import itertools

# Clear Climate Code
import read_config
from giss_data import valid, invalid, MISSING, BASE_YEAR
import parameters
import series


comb_log = open('log/comb.log','w')
pieces_log = open('log/pieces.log','w')

def comb_records(stream):
    """Combine records for the same station (the same id11) where
    possible.  For each station, the number of combined records will be
    at most the number of original records (in the case where no
    combining is possible), and each combined record is yielded."""

    return do_combine(stream, comb_log, get_best, combine)

def comb_pieces(stream):
    """comb_pieces() attempts to further combine the records produced
    by comb_records() - which have shorter overlaps - by comparing the
    annual anomalies of the years in which they do overlap, and
    finding ones for which the temperatures (in years which they do
    have in common) are on average closer together than the standard
    deviation of the combined record."""

    return do_combine(stream, pieces_log, get_longest, pieces_combine)

def do_combine(stream, log, select_func, combine_func):
    """Drive record combination.

    This is a filter driver function used by ``comb_records`` and
    ``comb_pieces``.

    :Param stream:
        The stream of records to filter.
    :Param log:
        Open log file file.
    :Param select_func:
        A function to call to select the 'best' record from a collection
        of records (belonging to the same station).
    :Param combine_func:
        A function to call to perform the data combining.

    """
    for id11, record_set in itertools.groupby(stream, lambda r: r.station_uid):
        log.write('%s\n' % id11)
        records = {}
        for record in record_set:
            records[record.uid] = record
            ann_mean, ann_anoms = series.monthly_annual(record.series)
            record.set_ann_anoms(ann_anoms)
            record.ann_mean = ann_mean
        ids = records.keys()
        while 1:
            if len(ids) == 1:
                yield records[ids[0]]
                break
            record_dict, begin, end = make_record_dict(records, ids)
            years = end - begin + 1
            record = select_func(record_dict)
            del record_dict[record.uid]
            sums, wgts = fresh_arrays(record, begin, years)
            log.write("\t%s %s %s -- %s\n" % (record.uid,
                record.first_valid_year(), record.last_valid_year(),
                record.source))
            combine_func(sums, wgts, begin, record_dict, log, record.uid)
            final_data = average(sums, wgts)
            record.set_series(begin * 12 + 1, final_data)
            yield record
            ids = record_dict.keys()
            if not ids:
                break

def combine(sums, wgts, begin, records, log, new_id_):
    while records:
        record, diff, overlap = get_longest_overlap(average(sums, wgts),
                                                    begin, records)
        if overlap < parameters.station_combine_min_overlap:
            log.write("\tno other records okay\n")
            return
        del records[record.uid]
        add(sums, wgts, diff, begin, record)
        log.write("\t %s %d %d %f\n" % (record.uid,
            record.first_valid_year(),
            record.last_valid_year(), diff))

def get_best(records):
    """Given a dict of records, return the "best" one.
    "best" considers the source of the record, preferring MCDW over
    USHCN over SUMOFDAY over UNKNOWN.

    (this is passed to do_combine() as a select_func argument)
    """

    ranks = {'MCDW': 4, 'USHCN2': 3, 'SUMOFDAY': 2, 'UNKNOWN': 1}
    best = 1
    longest = 0
    for record in sorted(records.values(), key=lambda r: r.uid):
        length = record.ann_anoms_good_count()
        rank = ranks[record.source]
        if rank > best:
            best = rank
            best_rec = record
        elif length > longest:
            longest = length
            longest_rec = record
    if best > 1:
        return best_rec
    return longest_rec

def pieces_combine(sums, wgts, begin, records, log, new_id):
    """The combine_func (passed to do_combine()) for comb_pieces().

    Combines remaining records that have insufficient overlap.
    """

    while records:
        record, diff_, overlap_ = get_longest_overlap(average(sums, wgts),
                                                      begin, records)
        log.write("\t %s %d %d\n" % (record.uid,
           record.first_valid_year(),
           record.last_valid_year()))

        is_okay = find_quintuples(sums, wgts, begin,
                                  record,
                                  new_id, log)

        if is_okay:
            del records[record.uid]
            add(sums, wgts, 0.0, begin, record)
        else:
            log.write("\t***no other pieces okay***\n")
            return

def get_longest(records):
    """Considering the records in the *records* dict, return the longest
    one.  This is the select_func (passed to do_combine()) used by
    comb_pieces."""

    def length(rec):
        """Length of a record, according to the number of valid annual
        anomalies."""

        return rec.ann_anoms_good_count()

    # :todo: If two records have the same "length", then which one is
    # chosen depends on the implementation of max() and
    # records.values(); it could matter.

    return max(records.values(), key=length)

def find_quintuples(new_sums, new_wgts, begin,
                    record, new_id, log):
    """Returns a boolean."""

    assert record.first_year == begin

    # An identifier common to all the log output.
    logid = "%s %s" % (new_id, record.uid)

    rec_begin = record.first_valid_year()
    rec_end = record.last_valid_year()

    actual_begin, actual_end = get_actual_endpoints(new_wgts, begin)

    max_begin = max(actual_begin, rec_begin)
    min_end = min(actual_end, rec_end)
    # Since max_begin and min_end are integers, this rounds fractional
    # middle years up.
    middle_year = int(.5 * (max_begin + min_end) + 0.5)
    log.write("max begin: %s\tmin end: %s\n" % (max_begin, min_end))

    new_data = average(new_sums, new_wgts)
    new_ann_mean, new_ann_anoms = series.monthly_annual(new_data)
    ann_std_dev = sigma(new_ann_anoms)
    log.write("ann_std_dev = %s\n" % ann_std_dev)
    offset = (middle_year - begin)
    new_len = len(new_ann_anoms)

    rec_ann_anoms = record.ann_anoms
    rec_ann_mean = record.ann_mean
    rec_len = len(rec_ann_anoms)

    # Whether we have an "overlap" or not.  We have an "overlap" if
    # within *rad* years either side of *middle_year* both records have
    # *parameters.station_combine_min_mid_year* valid annnual anomalies.
    ov_success = False
    # The overlap is "okay" when the difference in annual temperature is
    # below a certain threshold.
    okay_flag = False
    for rad in range(1, parameters.station_combine_bucket_radius + 1):
        count1 = sum1 = 0
        count2 = sum2 = 0
        # Try from -rad to rad (inclusive) around the middle year.
        for i in range(-rad, rad + 1):
            index = offset+i
            if 0 <= index < new_len:
                anom1 = new_ann_anoms[index]
                if valid(anom1):
                    sum1 += anom1 + new_ann_mean
                    count1 += 1
            if 0 <= index < rec_len:
                anom2 = rec_ann_anoms[index]
                if valid(anom2):
                    sum2 += anom2 + rec_ann_mean
                    count2 += 1
        if (count1 >= parameters.station_combine_min_mid_years
            and count2 >= parameters.station_combine_min_mid_years):
            log.write("overlap success: %s\n" % logid)
            ov_success = True
            avg1 = sum1 / float(count1)
            avg2 = sum2 / float(count2)
            diff = abs(avg1 - avg2)
            log.write("diff = %s\n" % diff)
            if diff < ann_std_dev:
                okay_flag = True
                log.write("combination success: %s\n" % logid)
            else:
                log.write("combination failure: %s\n" % logid)
            break
    if not ov_success:
        log.write("overlap failure: %s\n" % logid)
    log.write("counts: %d %d\n" % (count1, count2))
    return okay_flag


def adjust_helena(stream):
    """Modifies records as specified in config/combine_pieces_helena.in,
    by adding the delta to every datum for that station prior to the
    specified month.
    """
    helena_ds = read_config.get_helena_dict()
    for record in stream:
        id = record.uid
        if helena_ds.has_key(id):
            series = record.series
            this_year, month, summand = helena_ds[id]
            begin = record.first_year
            # Index of month specified by helena_ds
            M = (this_year - begin)*12 + month
            # All valid data up to and including M get adjusted
            for i in range(M+1):
                datum = series[i]
                if invalid(datum):
                    continue
                series[i] += summand
            record.set_series(record.first_month, series)
            del helena_ds[id]
        yield record

def drop_strange(data):
    """Drops station records, or parts of records, under control of
    the file 'config/Ts.strange.RSU.list.IN' file.
    """

    changes_dict = read_config.get_changes_dict()
    for record in data:
        changes = changes_dict.get(record.uid, [])
        series = record.series
        begin = record.first_year
        # :todo: Use record.last_year
        end = begin + (len(series)//12) - 1
        for (kind, year, x) in changes:
            if kind == 'years':
                # omit all the data from year1 to year2, inclusive
                year1 = year
                year2 = x
                if year1 <= begin and year2 >= end:
                    # Drop this whole record.  Note: avoids "else:"
                    # clause at end of "for" loop.
                    break

                # Clamp range of deleted years to the range of the
                # series.
                year1 = max(year1, begin)
                year2 = min(year2, end)
                if year2 < year1:
                    # Happens when deleted range is entirely outside the
                    # range of the series.  In which case, pass record
                    # unchanged.
                    continue
                # Invalidate the data.
                nmonths = (year2 + 1 - year1) * 12
                series[(year1-begin)*12:(year2+1-begin)*12] = [
                        MISSING] * nmonths

            else: # remove a single month
                series[(year-begin)*12 + x-1] = MISSING
        else:
            record.set_series(begin * 12 + 1, series)
            yield record


def alter_discont(data):
    """Modifies records as specified in config/Ts.discont.RS.alter.IN,
    by adding the delta to every datum for that station prior to the
    specified month.  Yes, this is very similar to adjust_helena().
    """

    alter_dict = read_config.get_alter_dict()
    for record in data:
        if alter_dict.has_key(record.uid):
            series = record.series
            (a_month, a_year, a_num) = alter_dict[record.uid]
            begin = record.first_year
            # Month index of the month in the config file.
            M = (a_year - begin)*12 + a_month - 1
            # Every (valid) month up to and not including the month in
            # question is adjusted.
            for i in range(M):
                if valid(series[i]):
                    series[i] += a_num
            record.set_series(record.first_month, series)

        yield record


def average(sums, counts):
    """Divide *sums* by *counts* to make a series of averages.
    Return an array with sums[i]/counts[i], and MISSING where
    counts[i] is zero.
    """

    assert len(sums) == len(counts)

    data = [MISSING] * (len(sums))

    for i,(sum,count) in enumerate(zip(sums, counts)):
        if count:
            data[i] = float(sum) / count

    return data

def add(sums, wgts, diff, begin, record):
    """Add the data from *record* to the *sums* and *wgts* arrays, first
    shifting it by subtracting *diff*.  *begin* is the first year for
    the *sums* and *wgts* array which is assumed (and asserted) to be
    the same as *record.first_year*."""

    assert record.first_year == begin
    for i,datum in enumerate(record.series):
        if invalid(datum):
            continue
        sums[i] += datum - diff
        wgts[i] += 1

def get_longest_overlap(target, begin, records):
    """Find the record in the *records* dict that has the longest
    overlap with the *target* by considering annual anomalies.  *target*
    is a sequence of monthly values starting in the year *begin*.

    A triple (record, diff, overlap) is returned; *diff* is the average
    difference in annual anomalies between *record* and *target*
    (positive when *record* is higher); *overlap* is the number of years
    in the overlap.  Even when there is no overlap _some_ record is
    returned and in that case *diff* is None and *overlap* is 0.
    
    Like other functions, assumes (and asserts) that *begin* is
    the first year for all the records.
    """

    # Annual mean, and annual anomaly sequence.
    mean, anoms = series.monthly_annual(target)
    overlap = 0
    diff = None
    # :todo: the records are consulted in an essentially arbitrary
    # order (chosen by the implementation of items()), but the order
    # may affect the result.  Tie breaks go to the last record consulted.
    for rec_id, record in records.items():
        assert record.first_year == begin
        common = [(rec_anom,anom)
          for rec_anom, anom in zip(record.ann_anoms, anoms)
          if valid(rec_anom) and valid(anom)]
        if len(common) < overlap:
            continue
        overlap = len(common)
        best_record = record
        S = sum((record.ann_mean+rec_anom) - (mean+anom)
                for rec_anom, anom in common)
        if common:
            diff = S / len(common)
    return best_record, diff, overlap

def make_record_dict(records, ids):
    """Build and return a fresh dictionary for a set of records.
    For each of the keys in *ids*, the corresponding entry in the
    *records* dictionary is consulted, and a new dictionary is made.

    (*record_dict*, *year_min*, *year_max*) is returned, where
    *record_dict* contains all the new dictionaries made; *year_min* and
    *year_max* are the minimum and maximum years with data, across all
    the records consulted.

    This function asserts that all the records have the same first year
    (which will be *year_min*)
    """

    record_dict = {}
    record_dict = dict((rec_id, records[rec_id]) for rec_id in ids)
    first_years = set(record.first_year for record in record_dict.values())
    assert 1 == len(first_years)
    y_min = list(first_years)[0]
    y_max = max(record.last_year for record in record_dict.values())
    return record_dict, y_min, y_max

def fresh_arrays(record, begin, years):
    """Make and return a fresh set of sums, and wgts arrays.  Each
    array is list (of length 12 * years).

    *begin* should be the starting year for the arrays, which is
    assumed, and asserted, to be the same as *record.first_year*.
    """

    assert record.first_year == begin

    nmonths = years * 12

    # Number of months in record.
    rec_months = len(record)
    assert rec_months <= nmonths

    sums = [0.0] * nmonths
    # Copy valid data rec_data into sums, assigning 0 for invalid data.
    sums[:rec_months] = (valid(x)*x for x in record.series)
    # Let wgts[i] be 1 where sums[i] is valid.
    wgts = [0] * nmonths
    wgts[:rec_months] = (int(valid(x)) for x in record.series)

    return sums, wgts


def sigma(list):
    # Remove invalid (missing) data.
    list = filter(valid, list)
    if len(list) == 0:
        return MISSING
    # Two pass method ensures argument to sqrt is always positive.
    mean = sum(list) / len(list)
    sigma_squared = sum((x-mean)**2 for x in list)
    return math.sqrt(sigma_squared/len(list))

def get_actual_endpoints(wgts, begin):
    """For the array of weights in *wgts* return the first and last
    calendar years that have some weight (contain a month with non-zero
    weight); assuming the array starts in year *begin*."""

    # Exact number of years.
    assert len(wgts) % 12 == 0
    y_min = 9999
    y_max = 0
    for i in range(0, len(wgts), 12):
        if sum(wgts[i:i+12]) > 0:
            y = i//12
            y_min = min(y_min, y)
            y_max = max(y_max, y)
    return begin+y_min, begin+y_max

def step1(record_source):
    """An iterator for step 1.  Produces a stream of
    `giss_data.Series` instances.

    :Param record_source:
        An iterable source of `giss_data.Series` instances (which it
        will assume are station records).

    """
    records = comb_records(record_source)
    helena_adjusted = adjust_helena(records)
    combined_pieces = comb_pieces(helena_adjusted)
    without_strange = drop_strange(combined_pieces)
    for record in alter_discont(without_strange):
        assert record.first_year == BASE_YEAR
        yield record

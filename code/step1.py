#! /usr/bin/env python
# $URL$
# $Rev$
#
# step1.py
#
# Nick Barnes, Ravenbrook Limited, 2008-08-06

"""
Python code reproducing the STEP1 part of the GISTEMP algorithm.

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
import struct
import itertools

import read_config
from giss_data import valid, invalid, MISSING
import parameters


def valid_mean(seq, min=1):
    """Takes a sequence, *seq*, and computes the mean of the valid
    items (using the valid() function).  If there are fewer than *min*
    valid items, the mean is MISSING."""

    count = 0
    sum = 0.0
    for x in seq:
        if valid(x):
            sum += x
            count += 1
    if count >= min:
        return sum/float(count)
    else:
        return MISSING

def monthly_anomalies(data):
    """Calculate monthly anomalies, by subtracting from every datum
    the mean for its month.  A pair of (monthly_mean, monthly_anom) is
    returned.  *monthly_mean* is a 12-long sequence giving the mean for
    each of the 12 months; *monthly_anom* is a 12-long sequence giving
    the anomalized series for each of the 12 months.
    
    The input data is a flat sequence, one datum per month.
    Effectively the data changes shape as it passes through this
    function.
    """

    years = len(data) // 12
    monthly_mean = []
    monthly_anom = []
    for m in range(12):
        row = data[m::12]
        mean = valid_mean(row)
        monthly_mean.append(mean)
        if valid(mean):
            def asanom(datum):
                """Convert a single datum to anomaly."""
                if valid(datum):
                    return datum - mean
                return MISSING
            monthly_anom.append(map(asanom, row))
        else:
            monthly_anom.append([MISSING]*years)
    return monthly_mean, monthly_anom

def monthly_annual(data):
    """Computing a mean and set of annual anomalies.  This has a
    particular algorithm (via monthly and then seasonal means and
    anomalies), which must be maintained for bit-for-bit compatibility
    with GISTEMP; maybe we can drop it later.
    """

    years = len(data) // 12
    monthly_mean, monthly_anom = monthly_anomalies(data)

    # :todo:
    # The seasonal calculation shoud be abstracted into a function.
    # inputs: years, monthly_anom, monthly_mean
    # outputs: seasonal_anom, seasonal_mean
    # Average monthly means to make seasonal means,
    # and monthly anomalies to make seasonal anomalies.
    seasonal_mean = []
    seasonal_anom = []
    # compute seasonal anomalies
    for s in range(4):
        if s == 0:
            months = [11, 0, 1]
        else:
            months = range(s*3-1, s*3+2)
        # need at least two valid months for a valid season
        seasonal_mean.append(valid_mean((monthly_mean[m] for m in months),
                                        min = 2))
        # A list of 3 data series, each being an extract for a
        # particular month
        month_in_season = []
        for m in months:
            row = monthly_anom[m] # the list of anomalies for month m
            if m == 11:
                # For December, we take the December of the previous
                # year.  Which we do by offsetting the array, and not
                # using the most recent December.
                row[1:] = row[:-1]
                row[0] = MISSING
            month_in_season.append(row)
        seasonal_anom_row = []
        for n in range(years):
            m = valid_mean((month_in_season[i][n] for i in range(3)),
                           min = 2)
            seasonal_anom_row.append(m)
        seasonal_anom.append(seasonal_anom_row)

    # Average seasonal means to make annual mean,
    # and average seasonal anomalies to make annual anomalies
    # (note: annual anomalies are December-to-November).
    annual_mean = valid_mean(seasonal_mean, min = 3)
    annual_anom = []
    for n in range(years):
        annual_anom.append(valid_mean((seasonal_anom[s][n] for s in range(4)),
                                      min = 3))
    return (annual_mean, annual_anom)


def average(sums, counts, years):
    """Return an array with sums[i]/counts[i], and MISSING where
    counts[i] is zero.
    """

    assert len(sums) == years * 12
    assert len(counts) == years * 12

    data = [MISSING] * (years*12)

    for i in range(len(counts)):
        count = counts[i]
        if count != 0:
            data[i] = float(sums[i]) / count

    return data

def final_average(sums, wgts, years, begin):
    """Compute an average from *sums* and *wgts*, return a new
    (*begin*,*data*) pair.  The results are trimmed so that the first
    and last years have valid data.
    """

    y_min, y_max = 9999, -9999
    assert len(wgts) == 12*years
    # Find the minimum and maximum years with data.
    for i in range(years):
        if sum(wgts[i*12:(i+1)*12]) > 0:
            y_min = min(y_min, i)
            y_max = max(y_max, i)
    if y_min == 0 and y_max == years - 1:
        return begin, average(sums, wgts, years)
    years = y_max - y_min + 1
    begin = begin + y_min
    month_base = y_min * 12
    month_limit = (y_max+1) * 12
    sums = sums[month_base:month_limit]
    wgts = wgts[month_base:month_limit]
    assert len(wgts) == 12*years
    return begin, average(sums, wgts, years)

def add(sums, wgts, diff, begin, record):
    """Add the data from *record* to the *sums* and *wgts* arrays, first
    shifting it by subtracting *diff*."""

    rec_begin = record.first_year
    rec_years = record.last_year - record.first_year + 1
    rec_data = record.series
    assert len(rec_data) == 12*rec_years
    offset = rec_begin - begin
    offset *= 12
    for i in range(len(rec_data)):
        datum = rec_data[i]
        if invalid(datum):
            continue
        index = i + offset
        sums[index] += datum - diff
        wgts[index] += 1

def get_longest_overlap(new_data, begin, records):
    """Find the record in the *records* dict that has the longest
    overlap with the *new_data* by considering annual anomalies.
    """

    ann_mean, ann_anoms = monthly_annual(new_data)
    overlap = 0
    # :todo: the records are consulted in an essentially arbitrary
    # order (chosen by the implementation of items()), but the order
    # may affect the result.
    # Tie breaks go to the last record consulted.
    for rec_id, record in records.items():
        rec_ann_anoms = record.ann_anoms
        rec_ann_mean = record.ann_mean
        rec_years = record.last_year - record.first_year + 1
        rec_begin = record.first_year
        sum = wgt = 0
        for n in range(rec_years):
            rec_anom = rec_ann_anoms[n]
            if invalid(rec_anom):
                continue
            year = n + rec_begin
            anom = ann_anoms[year - begin]
            if invalid(anom):
                continue
            wgt += 1
            sum += (rec_ann_mean + rec_anom) - (ann_mean + anom)
        if wgt < parameters.station_combine_min_overlap:
            continue
        if wgt < overlap:
            continue
        overlap = wgt
        diff = sum / wgt
        best_id = rec_id
        best_record = record
    if overlap < parameters.station_combine_min_overlap:
        return 0, 0, MISSING
    return best_record, best_id, diff

def combine(sums, wgts, begin, years, records, log, new_id=None):
    while records:
        record, rec_id, diff = get_longest_overlap(average(sums, wgts, years),
                                                   begin, records)
        if invalid(diff):
            log.write("\tno other records okay\n")
            return
        del records[rec_id]
        add(sums, wgts, diff, begin, record)
        log.write("\t %s %d %d %f\n" % (rec_id, record.first_year,
                                        record.last_year - 1, diff))

def get_best(records):
    """Given a set of records (a dict really), return the best one, and
    its key in the *records* dict.
    """

    ranks = {'MCDW': 4, 'USHCN2': 3, 'SUMOFDAY': 2, 'UNKNOWN': 1}
    best = 1
    longest = 0
    for rec_id in sorted(records.keys()):
        record = records[rec_id]
        length = record.ann_anoms_good_count
        rank = ranks[record.source]
        if rank > best:
            best = rank
            best_rec = record
            best_id = rec_id
        elif length > longest:
            longest = length
            longest_rec = record
            longest_id = rec_id
    if best > 1:
        return best_rec, best_id
    return longest_rec, longest_id

def make_record_dict(records, ids):
    """Build and return a fresh dictionary for a set of records.
    For each of the keys in *ids*, the corresponding entry in the
    *records* dictionary is consulted, and a new dictionary is made.

    (*record_dict*, *year_min*, *year_max*) is returned, where
    *record_dict* contains all the new dictionaries made; *year_min* and
    *year_max* are the minimum and maximum years with data, across all
    the records consulted.
    """

    record_dict = {}
    y_min, y_max = 9999, -9999
    for rec_id in ids:
        record = records[rec_id]
        begin = record.first_year
        end = record.last_year
        y_min = min(y_min, begin)
        y_max = max(y_max, end)
        record_dict[rec_id] = record
    return record_dict, y_min, y_max

def fresh_arrays(record, begin, years):
    """Make and return a fresh set of sums, and wgts arrays.  Each
    array is list (of length 12 * years).

    *begin* should be the starting year for the arrays, which must
    be no later than the starting year for the record.
    """

    nmonths = years * 12

    rec_data = record.series
    rec_begin = record.first_year
    rec_years = record.last_year - record.first_year + 1
    # Number of months in record.
    rec_months = rec_years * 12
    assert rec_months == record.n
    assert rec_months == len(rec_data)
    # The record may begin at a later year from the arrays we are
    # creating, so we need to offset it when we copy.
    offset = rec_begin - begin
    assert offset >= 0
    offset *= 12

    sums = [0.0] * nmonths
    # Copy valid data rec_data into sums, assigning 0 for invalid data.
    sums[offset:offset+rec_months] = (valid(x)*x for x in rec_data)
    # Let wgts[i] be 1 where sums[i] is valid.
    wgts = [0] * nmonths
    wgts[offset:offset+rec_months] = (int(valid(x)) for x in rec_data)

    return sums, wgts


def comb_records(stream):
    """Combine records for the same station (the same id11) where
    possible.  For each station, the number of combined records will be
    at most the number of original records (in the case where no
    combining is possible), and each combined record is yielded."""

    return do_combine(stream, open('log/comb.log','w'), get_best, combine)


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

def sigma(list):
    # Remove invalid (missing) data.
    list = filter(valid, list)
    if len(list) == 0:
        return MISSING
    # Two pass method ensures argument to sqrt is always positive.
    mean = sum(list) / len(list)
    sigma_squared = sum((x-mean)**2 for x in list)
    return math.sqrt(sigma_squared/len(list))

# Annoyingly similar to get_longest_overlap
def pieces_get_longest_overlap(new_data, begin, records):
    ann_mean, ann_anoms = monthly_annual(new_data)
    overlap = 0
    for rec_id, record in records.items():
        rec_ann_anoms = record.ann_anoms
        rec_years = record.last_year - record.first_year + 1
        rec_begin = record.first_year
        wgt = 0
        for n in range(rec_years):
            rec_anom = rec_ann_anoms[n]
            if invalid(rec_anom):
                continue
            year = n + rec_begin
            anom = ann_anoms[year - begin]
            if invalid(anom):
                continue
            wgt = wgt + 1
        if wgt < overlap:
            continue
        overlap = wgt
        best_id = rec_id
        best_record = record
    return best_record, best_id

# :todo: unify with code in final_average
def get_actual_endpoints(wgts, begin, years):
    assert len(wgts) == 12*years
    y_min = 9999
    y_max = 0
    for i in range(years):
        if sum(wgts[i*12:(i+1)*12]) > 0:
            y_min = min(y_min, i)
            y_max = max(y_max, i)
    return begin+y_min, begin+y_max

def find_quintuples(new_sums, new_wgts,
                    begin, years, record, rec_begin,
                    new_id, rec_id, log):
    rec_begin = record.first_year
    rec_end = rec_begin + record.last_year - record.first_year

    actual_begin, actual_end = get_actual_endpoints(new_wgts, begin, years)

    max_begin = max(actual_begin, rec_begin)
    min_end = min(actual_end, rec_end)
    middle_year = int(.5 * (max_begin + min_end) + 0.5)
    log.write("max begin: %s\tmin end: %s\n" % (max_begin, min_end))

    new_data = average(new_sums, new_wgts, years)
    new_ann_mean, new_ann_anoms = monthly_annual(new_data)
    ann_std_dev = sigma(new_ann_anoms)
    log.write("ann_std_dev = %s\n" % ann_std_dev)
    new_offset = (middle_year - begin)
    new_len = len(new_ann_anoms)

    rec_ann_anoms = record.ann_anoms
    rec_ann_mean = record.ann_mean
    rec_offset = (middle_year - rec_begin)
    rec_len = len(rec_ann_anoms)

    ov_success = 0
    okay_flag = 0
    for rad in range(1, parameters.station_combine_bucket_radius + 1):
        count1 = sum1 = 0
        count2 = sum2 = 0
        for i in range(0, rad + 1):
            for sign in [-1, 1]:
                if sign == 1 and i == 0:
                    continue
                index1 = i * sign + new_offset
                index2 = i * sign + rec_offset
                if index1 < 0 or index1 >= new_len:
                    anom1 = MISSING
                else:
                    anom1 = new_ann_anoms[index1]
                if index2 < 0 or index2 >= rec_len:
                    anom2 = MISSING
                else:
                    anom2 = rec_ann_anoms[index2]
                if valid(anom1):
                    sum1 += anom1 + new_ann_mean
                    count1 += 1
                if valid(anom2):
                    sum2 += anom2 + rec_ann_mean
                    count2 += 1
        if (count1 >= parameters.station_combine_min_mid_years
            and count2 >= parameters.station_combine_min_mid_years):
            log.write("overlap success: %s %s\n" % (new_id, rec_id))
            ov_success = 1
            avg1 = sum1 / float(count1)
            avg2 = sum2 / float(count2)
            diff = abs(avg1 - avg2)
            log.write("diff = %s\n" % diff)
            if diff < ann_std_dev:
                okay_flag = 1
                log.write("combination success: %s %s\n" % (new_id, rec_id))
            else:
                log.write("combination failure: %s %s\n" % (new_id, rec_id))
            break
    if not ov_success:
        log.write("overlap failure: %s %s\n" % (new_id, rec_id))
    log.write("counts: %s\n" % ((count1, count2),))
    return okay_flag

def pieces_combine(sums, wgts, begin, years, records, log, new_id):
    while records:
        record, rec_id = pieces_get_longest_overlap(average(sums, wgts, years),
                                                    begin, records)
        rec_begin = record.first_year
        rec_end = rec_begin + record.last_year - record.first_year

        log.write("\t %s %d %d\n" % (rec_id, rec_begin, rec_end - 1))

        is_okay = find_quintuples(sums, wgts,
                                  begin, years, record, rec_begin,
                                  new_id, rec_id, log)

        if is_okay:
            del records[rec_id]
            add(sums, wgts, 0.0, begin, record)
            log.write("\t %s %d %d\n" % (rec_id, rec_begin, rec_end - 1))
        else:
            log.write("\t***no other pieces okay***\n")
            return

def get_longest(records):
    """Considering the records in the *records* dict, return the longest
    one."""

    longest = 0
    # :todo: The order depends on the implementation of records.items,
    # and could matter.
    for rec_id, record in records.items():
        length = record.ann_anoms_good_count
        if length > longest:
            longest = length
            longest_rec = record
            longest_id = rec_id
    return longest_rec, longest_id

def do_combine(stream, log, select_func, combine_func):
    """Drive record combination.

    This is a filter driver function used by ``comb_records`` and
    ``comb_pieces``.

    :Param stream:
        The stream of records to filter.
    :Param log:
        Open log file file.
    :Param select_func:
        A function to call to select the the 'best' records from a set of
        records belonging to the same station.
    :Param combine_func:
        A function to call to perform the data combining.

    """
    for id11, record_set in itertools.groupby(stream, lambda r: r.station_uid):
        log.write('%s\n' % id11)
        records = {}
        for record in record_set:
            records[record.uid] = record
            ann_mean, ann_anoms = monthly_annual(record.series)
            record.set_ann_anoms(ann_anoms)
            record.ann_mean = ann_mean
        ids = records.keys()
        while 1:
            if len(ids) == 1:
                yield records[ids[0]]
                break
            record_dict, begin, end = make_record_dict(records, ids)
            years = end - begin + 1
            record, rec_id = select_func(record_dict)
            del record_dict[rec_id]
            sums, wgts = fresh_arrays(record, begin, years)
            log.write("\t%s %s %s -- %s\n" % (rec_id, begin,
                begin + years -1, record.source))
            combine_func(sums, wgts, begin, years, record_dict, log, rec_id)
            begin, final_data = final_average(sums, wgts, years, begin)
            record.set_series(begin * 12 + 1, final_data)
            yield record
            ids = record_dict.keys()
            if not ids:
                break


def comb_pieces(stream):
    """comb_pieces() attempts to further combine the records produced
    by comb_records() - which have shorter overlaps - by comparing the
    annual anomalies of the years in which they do overlap, and
    finding ones for which the temperatures (in years which they do
    have in common) are on average closer together than the standard
    deviation of the combined record."""

    return do_combine(stream, open('log/pieces.log','w'),
                      get_longest, pieces_combine)


def drop_strange(data):
    """Drops station records, or parts of records, under control of
    the file 'config/Ts.strange.RSU.list.IN' file.
    """

    changes_dict = read_config.get_changes_dict()
    for record in data:
        changes = changes_dict.get(record.uid, [])
        series = record.series
        begin = record.first_year
        end = begin + (len(series)//12) - 1
        for (kind, year, x) in changes:
            if kind == 'years':
                # omit all the data from year1 to year2, inclusive
                year1 = year
                year2 = x
                if year1 <= begin and year2 >= end:
                    # drop this whole record
                    break

                if year1 <= begin:
                    if year2 >= begin:
                        # trim at the start
                        series = series[(year2 - begin + 1)*12:]
                        begin = year2 + 1
                    continue

                if year2 >= end:
                    if year1 <= end:
                        # trim at the end
                        series = series[:(year1 - begin)*12]
                        end = year1 - 1
                    continue
                # remove some years from mid-series
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


def step1(record_source):
    """An iterator for step 1.  Produces a stream of
    `giss_data.StationRecord` instances.

    :Param record_source:
        An iterable source of `giss_data.StationRecord` instances.

    """
    records = comb_records(record_source)
    helena_adjusted = adjust_helena(records)
    combined_pieces = comb_pieces(helena_adjusted)
    without_strange = drop_strange(combined_pieces)
    for record in alter_discont(without_strange):
        yield record

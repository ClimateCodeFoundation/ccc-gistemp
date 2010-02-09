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
import array
import itertools

import read_config
import giss_data


# TODO: Make valid/invalid checking a giss_data function.
def invalid(x):
    """Test for invalid datum ("equal" to the BAD value, for some
    definition of "equal").
    """

    # If you're feeling spooky about the BAD value, re-enable this:
    # if abs(x-BAD) < 0.1:
    #     assert x == BAD

    return x == giss_data.XMISSING

def valid(x):
    """Test for valid datum.  See invalid()."""

    # It's important that this obey Iverson's convention: in other words
    # return 0 or 1 (or False or True, which it does).
    return not invalid(x)

def sum_valid(seq):
    """Takes a sequence, *seq*, and sums up all the valid items (using
    the valid() function).  *sum*,*count* is returned (where *count* is
    the number of valid items).
    """

    count = 0
    a = 0.0
    for x in seq:
        if valid(x):
            a += x
            count += 1
    return a,count


def month_anomaly(data):
    """Convert data to monthly anomalies, by subtracting from every
    datum the mean for its month.  A pair of (monthly_mean,monthly_anom)
    is returned.
    """

    years = len(data) // 12
    #print "monthly_annual-1", years
    # For each of the 12 months, its mean.
    monthly_mean = [giss_data.MISSING]*12
    # Every datum converted to anomaly by subtracting the mean
    # for its month.
    monthly_anom = array.array('d', data)
    for m in range(12):
        # *row* contains the data for one month of the year.
        row = data[m::12]
        sum,count = sum_valid(row)
        #print "monthly_annual-2", sum,count, row
        if count > 0:
            monthly_mean[m] = sum / float(count)
        mean = monthly_mean[m]
        def asanom(datum):
            """Convert a single datum to anomaly."""
            if valid(datum):
                return datum - mean
            return giss_data.MISSING
        if valid(mean):
            anom_row = array.array('d', map(asanom, row))
        else:
            anom_row = array.array('d',[giss_data.MISSING]*years)
        monthly_anom[m::12] = anom_row
    return monthly_mean, monthly_anom

def monthly_annual(data):
    """Computing a mean and set of annual anomalies.  This has a
    particular algorithm (via monthly and then seasonal means and
    anomalies), which must be maintained for bit-for-bit compatibility
    with GISTEMP; maybe we can drop it later.
    """

    years = len(data) // 12
    monthly_mean, monthly_anom = month_anomaly(data)

    # :todo:
    # The seasonal calculation shoud be abstracted into a function.
    # inputs: years, monthly_anom, monthly_mean
    # outputs: seasonal_anom, seasonal_mean
    # Average monthly means to make seasonal means,
    # and monthly anomalies to make seasonal anomalies.
    seasonal_mean = array.array('d', [giss_data.MISSING] * 4)
    seasonal_anom = []
    # compute seasonal anomalies
    for s in range(4):
        if s == 0:
            months = [11, 0, 1]
        else:
            months = range(s*3-1, s*3+2)
        sum,count = sum_valid(monthly_mean[m] for m in months)
        # need at least two valid months for a valid season
        if count > 1:
            seasonal_mean[s] = sum / float(count)
        seasonal_anom_row = array.array('d', [giss_data.MISSING] * years)
        # A list of 3 data series, each being an extract for a
        # particular month
        month_in_season = []
        for m in months:
            row = monthly_anom[m::12]
            if m == 11:
                # For december, we take the december of the previous
                # year.  Which we do by offsetting the array, and not
                # using the most recent december.
                row[1:] = row[:-1]
                row[0] = giss_data.MISSING
            month_in_season.append(row)
        for n in range(years):
            sum = 0.0
            count = 0
            for i in range(3):
                m_anom = month_in_season[i][n]
                if valid(m_anom):
                    sum += m_anom
                    count += 1
            # need at least two valid months for a valid season
            if count > 1:
                seasonal_anom_row[n] = sum / float(count)
        seasonal_anom.append(seasonal_anom_row)
    
    # Average seasonal means to make annual mean,
    # and average seasonal anomalies to make annual anomalies
    # (note: annual anomalies are December-to-November).
    sum,count = sum_valid(seasonal_mean)
    # need 3 valid seasons for a valid year
    if count > 2:
        annual_mean = sum / float(count)
    else:
        annual_mean = giss_data.MISSING
    annual_anom = array.array('d',[giss_data.MISSING]*years)
    for n in range(years):
        sum,count = sum_valid(seasonal_anom[s][n] for s in range(4))
        # need 3 valid seasons for a valid year
        if count > 2:
            annual_anom[n] = sum / float(count)
    return (annual_mean, annual_anom)


MIN_OVERLAP = 4
comb_log = None

def average(sums, counts, out_data, years):
    """Where the *counts* array indicates that *sums* has a valid sum,
    update the corresponding entry in *out_data* with the average
    (sum/count).

    This function is scheduled for termination, see fresh_average.
    """

    assert len(sums) == years * 12
    assert len(counts) == years * 12
    assert len(out_data) == years * 12

    for i in range(len(counts)):
        count = counts[i]
        if count != 0:
            out_data[i] = sums[i] / count

def fresh_average(sums, counts, years):
    """Same as average, but returns a fresh data array.  The plan is to
    eliminate average by replacing calls to it with calls to this
    function.  Then rename this as average."""

    data = [giss_data.MISSING] * (years*12)
    average(sums, counts, data, years)
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
        return begin, fresh_average(sums, wgts, years)
    years = y_max - y_min + 1
    begin = begin + y_min
    month_base = y_min * 12
    month_limit = (y_max+1) * 12
    sums = sums[month_base:month_limit]
    wgts = wgts[month_base:month_limit]
    assert len(wgts) == 12*years
    return begin, fresh_average(sums, wgts, years)

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

def get_longest_overlap(sums, wgts, begin, years, records):
    """Find the record in the *records* dict that has the longest
    overlap with the accumulated data in *sums* and *wgts* by
    considering annual anomalies.  An overlap has to be at least
    MIN_OVERLAP years to count.
    """

    #print "get_longest_overlap"
    new_data = fresh_average(sums, wgts, years)
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
        if wgt < MIN_OVERLAP:
            continue
        if wgt < overlap:
            continue
        overlap = wgt
        diff = sum / wgt
        best_id = rec_id
        best_record = record
    if overlap < MIN_OVERLAP:
        return 0, 0, giss_data.MISSING
    return best_record, best_id, diff

def combine(sums, wgts, begin, years, records):
    while records:
        record, rec_id, diff = get_longest_overlap(sums, wgts, begin,
                years, records)
        if invalid(diff):
            comb_log.write("\tno other records okay\n")
            return
        del records[rec_id]
        add(sums, wgts, diff, begin, record)
        rec_begin = record.first_year
        comb_log.write("\t %s %d %d %f\n" %
          (rec_id, rec_begin, record.last_year - record.first_year + rec_begin - 1, diff))

def get_best(records):
    """Given a set of records (a dict really), return the best one, and
    its key in the *records* dict.
    """

    ranks = {'MCDW': 4, 'USHCN2': 3, 'SUMOFDAY': 2, 'UNKNOWN': 1}
    best = 1
    longest = 0
    for rec_id in sorted(records.keys()):
        record = records[rec_id]
        source = record.source
        length = record.ann_anoms_good_count
        rank = ranks[source]
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
        ann_mean, ann_anoms = monthly_annual(record.series)
        # Let *length* be the number of valid data in ann_anoms.
        #print rec_id, ann_mean, len(ann_anoms), ':', ann_anoms[:5]
        length = 0
        for anom in ann_anoms:
            length += valid(anom)
        record_dict[rec_id] = giss_data.StationRecord(record.uid)
        record_dict[rec_id].set_ann_anoms(ann_anoms)
        record_dict[rec_id].ann_mean = ann_mean
        record_dict[rec_id].set_series_from_tenths(record.first_month,
            record.series_as_tenths)
        record_dict[rec_id].source = record.source
    return record_dict, y_min, y_max

def fresh_arrays(record, begin, years):
    """Make and return a fresh set of sums, and wgts arrays.  Each
    array is list (of length 12 * years).

    *begin* should be the starting year for the arrays, which must
    be no later than the starting year for the record.
    """

    nmonths = years * 12

    rec_data = record.series
    rec_begin, rec_years = record.first_year, record.last_year - record.first_year + 1
    # Number of months in record.
    rec_months = rec_years * 12
    assert rec_months == record.n
    assert rec_months == len(rec_data)
    # The record may begin at a later year from the arrays we are
    # creating, so we need to offset it when we copy.
    offset = rec_begin - begin
    assert offset >= 0
    offset *= 12

    try:
        sums = [0.0] * nmonths
    except MemoryError:
        print "???", nmonths
        raise
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

    global comb_log
    comb_log = open('log/comb.log','w')

    for id11, record_set in itertools.groupby(stream, lambda r: r.station_uid):
        comb_log.write('%s\n' % id11)
        records = {}
        for record in record_set:
            records[record.uid] = record
        ids = records.keys()
        while 1:
            if len(ids) == 1:
                yield records[ids[0]]
                break
            record_dict, begin, end = make_record_dict(records, ids)
            years = end - begin + 1
            record, rec_id = get_best(record_dict)
            del record_dict[rec_id]
            sums, wgts = fresh_arrays(record, begin, years)
            new_record = record.copy()
            new_record.source = record.source
            new_record.ann_mean = record.ann_mean
            new_record.ann_anoms = record.ann_anoms

            comb_log.write ("\t%s %s %s -- %s\n" % (rec_id, begin, end,
                record.source))
            combine(sums, wgts, begin, years, record_dict)
            begin,final_data = final_average(sums, wgts, years, begin)
            new_record.set_series(begin * 12 + 1, final_data)
            yield new_record
            ids = record_dict.keys()
            if not ids:
                break


def adjust_helena(stream):
    """Modifies records as specified in config/combine_pieces_helena.in,
    by adding the delta to every datum for that station prior to the
    specified month.  Returns an iterator.
    """
    helena_ds = read_config.get_helena_dict()
    for record in stream:
        id = record.uid
        if helena_ds.has_key(id):
            series = list(record.series)
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

MIN_MID_YEARS = 5            # MIN_MID_YEARS years closest to ymid
BUCKET_RADIUS = 10
pieces_log = None

def sigma(list):
    sos = sum = count = 0
    for x in list:
        if invalid(x):
            continue
        sos += x * x
        sum += x
        count = count + 1
    if count == 0:
        return giss_data.MISSING
    else:
        mean = sum/count
        return math.sqrt(sos / count - mean * mean)

# Annoyingly similar to get_longest_overlap
def pieces_get_longest_overlap(sums, wgts, begin, years, records):
    new_data = fresh_average(sums, wgts, years)
    ann_mean, ann_anoms = monthly_annual(new_data)
    overlap = 0
    length = 0
    for rec_id, record in records.items():
        rec_ann_anoms = record.ann_anoms
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
            wgt = wgt + 1
            sum = sum + rec_anom - anom
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
                    new_id, rec_id):
    rec_begin = record.first_year
    rec_end = rec_begin + record.last_year - record.first_year

    actual_begin, actual_end = get_actual_endpoints(new_wgts, begin, years)

    max_begin = max(actual_begin, rec_begin)
    min_end = min(actual_end, rec_end)
    middle_year = int(.5 * (max_begin + min_end) + 0.5)
    pieces_log.write("max begin: %s\tmin end: %s\n" % (max_begin, min_end))

    new_data = fresh_average(new_sums, new_wgts, years)
    new_ann_mean, new_ann_anoms = monthly_annual(new_data)
    ann_std_dev = sigma(new_ann_anoms)
    pieces_log.write("ann_std_dev = %s\n" % ann_std_dev)
    new_offset = (middle_year - begin)
    new_len = len(new_ann_anoms)

    rec_ann_anoms = record.ann_anoms
    rec_ann_mean = record.ann_mean
    rec_offset = (middle_year - rec_begin)
    rec_len = len(rec_ann_anoms)

    ov_success = 0
    okay_flag = 0
    for rad in range(1, BUCKET_RADIUS + 1):
        count1 = sum1 = 0
        count2 = sum2 = 0
        for i in range(0, rad + 1):
            for sign in [-1, 1]:
                if sign == 1 and i == 0:
                    continue
                index1 = i * sign + new_offset
                index2 = i * sign + rec_offset
                if index1 < 0 or index1 >= new_len:
                    anom1 = giss_data.MISSING
                else:
                    anom1 = new_ann_anoms[index1]
                if index2 < 0 or index2 >= rec_len:
                    anom2 = giss_data.MISSING
                else:
                    anom2 = rec_ann_anoms[index2]
                if valid(anom1):
                    sum1 += anom1 + new_ann_mean
                    count1 += 1
                if valid(anom2):
                    sum2 += anom2 + rec_ann_mean
                    count2 += 1
        if count1 >= MIN_MID_YEARS and count2 >= MIN_MID_YEARS:
            pieces_log.write("overlap success: %s %s\n" % (new_id, rec_id))
            ov_success = 1
            avg1 = sum1 / float(count1)
            avg2 = sum2 / float(count2)
            diff = abs(avg1 - avg2)
            pieces_log.write("diff = %s\n" % diff)
            if diff < ann_std_dev:
                okay_flag = 1
                pieces_log.write("combination success: %s %s\n" % (new_id, rec_id))
            else:
                pieces_log.write("combination failure: %s %s\n" % (new_id, rec_id))
            break
    if not ov_success:
        pieces_log.write("overlap failure: %s %s\n" % (new_id, rec_id))
    pieces_log.write("counts: %s\n" % ((count1, count2),))
    return okay_flag

def pieces_combine(sums, wgts, begin, years, records, new_id):
    while records:
        record, rec_id = pieces_get_longest_overlap(sums, wgts, begin, years, records)
        rec_begin = record.first_year
        rec_end = rec_begin + record.last_year - record.first_year

        pieces_log.write("\t %s %d %d\n" % (rec_id, rec_begin, record.last_year - record.first_year + rec_begin - 1))

        is_okay = find_quintuples(sums, wgts,
                                  begin, years, record, rec_begin,
                                  new_id, rec_id)

        if is_okay:
            del records[rec_id]
            add(sums, wgts, 0.0, begin, record)
            pieces_log.write("\t %s %d %d\n" % (rec_id, rec_begin, record.last_year - record.first_year + rec_begin - 1))
        else:
            pieces_log.write("\t***no other pieces okay***\n")
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

def comb_pieces(stream):
    global pieces_log
    pieces_log = open('log/pieces.log','w')

    # TODO: Very similar structure to comb_records
    for id11, record_set in itertools.groupby(stream, lambda r: r.station_uid):
        pieces_log.write('%s\n' % id11)
        records = {}
        for record in record_set:
            records[record.uid] = record
            stid = record.station_uid
        ids = records.keys()
        while 1:
            if len(ids) == 1:
                yield records[ids[0]]
                break
            record_dict, begin, end = make_record_dict(records, ids)
            record, rec_id = get_longest(record_dict)
            years = end - begin + 1
            del record_dict[rec_id]
            sums, wgts = fresh_arrays(record, begin, years)
            new_record = record.copy()
            new_record.source = record.source
            new_record.ann_mean = record.ann_mean
            new_record.ann_anoms = record.ann_anoms

            pieces_log.write("\t%s %s %s -- %s\n" % (rec_id, begin,
                begin + years -1, record.source))
            pieces_combine(sums, wgts, begin, years, record_dict, rec_id)
            begin, final_data = final_average(sums, wgts, years, begin)
            new_record.set_series(begin * 12 + 1, final_data)
            yield new_record
            ids = record_dict.keys()
            if not ids:
                break

def drop_strange(data):
    """Drops data from station records, under control of the file
    'config/Ts.strange.RSU.list.IN' file.  Returns an iterator.
    """

    changes_dict = read_config.get_changes_dict()
    for record in data:
        changes = changes_dict.get(record.uid, [])
        series = list(record.series)
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
                        #record.first_year = begin
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
                        giss_data.MISSING] * nmonths

            else: # remove a single month
                series[(year-begin)*12 + x-1] = giss_data.MISSING

        else:
            record.set_series(begin * 12 + 1, series)
            yield record


def alter_discont(data):
    """Modifies records as specified in config/Ts.discont.RS.alter.IN,
    by adding the delta to every datum for that station prior to the
    specified month.  Returns an iterator.
    """

    alter_dict = read_config.get_alter_dict()
    for record in data:
        if alter_dict.has_key(record.uid):
            series = list(record.series)
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


class Step1Iterator(object):
    """An iterator for step 1.
    
    An instance of this class acts as an iterator that produces a stream of
    `giss_data.StationRecord` instances.

    """
    def __init__(self, record_source):
        """Constructor:

        :Param record_source:
            An iterable source of `giss_data.StationRecord` instances.

        """
        self.record_source = record_source

    def __iter__(self):
        return self._it()

    def _it(self):
        records = comb_records(self.record_source)
        helena_adjusted = adjust_helena(records)
        combined_pieces = comb_pieces(helena_adjusted)
        without_strange = drop_strange(combined_pieces)
        for record in alter_discont(without_strange):
            yield record

    
def step1(record_source):
    return Step1Iterator(record_source)


# Notes:
#
# - comb_records() is a filter, combining several separate "scribal"
# station records into a single record, based on a minimum overlap
# (number of years with both stations valid) between pairs of records,
# attempting to make single coherent records.
#
# - adjust_helena() makes discontinuity adjustments to particular
# records, under control of a configuration file.  In fact, the
# configuration file specifies a single adjustment to a single record
# (St Helena).
#
# - comb_pieces() is another filter, which further attempts to combine
# the records produced by comb_records() - which therefore have
# shorter overlaps - by comparing the annual anomalies of the years in
# which they do overlap, and finding ones for which the temperatures
# (in years which they do have in common) are on average closer
# together than the standard deviation of the combined record.  I
# think.
# 
# - drop_strange() is another filter, which discards some station
# records, or parts of records, under control of a configuration file.
# 
# - alter_discont() makes discontinuity adjustments to records under
# control of a configuration file.  Yes, this is very similar to the
# St Helena behaviour of comb_pieces().
#
# - write_to_file() writes the results to work/Ts.txt

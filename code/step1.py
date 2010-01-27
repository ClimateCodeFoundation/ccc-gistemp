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

def round_series(series):
    """Round every element in *series*, in-place, to the nearest 0.1.
    Returns *series*.
    """

    for i in range(len(series)):
        series[i] = float(math.floor(series[i] * 10.0 + 0.5)) * 0.1
    return series

# BAD is a value used when the datum is not valid.  It is approximately
# 999.9 (9999 in tenths of degree C).  The BAD value only originates in
# data from two places.  The function from_lines inserts into the data
# where it detects data marked as invalid in the input (which uses a
# different invalid marker); the function round_series rounds all data
# (including any BAD value) to the nearest 0.1.  If we can ensure that
# the BAD value remains unchanged when passing through round_series then
# the invalid() test, below, can use an exact equality test.
# Defining BAD like this makes it insensitive to the exact definition of
# round_series.
BAD = round_series([999.9])[0]
assert round_series([BAD]) == [BAD]

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

# For each station record we carry around a dict and a series of
# temperature records.  The series is a linear array of floating-point
# temperatures: data[m], where m starts at 0 for the January of the
# first year in the record, and increments for every month.
# Missing values are recorded as 999.9 (ish).
#
# The dict has various fields, including the following:
# 'id'      : the station ID
# 'begin'   : first year of data
# 'name'
# 'lat'     : (floating-point) latitude
# 'lon'     : (floating-point) longitude
# 'elevs'
# 'elevg'
# 'pop'
# 'ipop'
# 'topo'
# 'stveg',
# 'stloc'
# 'iloc'
# 'airstn'
# 'itowndis'
# 'grveg'
# 
# The data and dict items are originally created (from the
# v2.mean_comb file and metadata files) by the from_lines() and
# read_v2() functions.

def from_lines(lines):
    """*lines* is a list of lines (strings) that comprise a station's
    entire record.  The lines are an extract from a file in the same
    format as the GHCN file v2.mean.Z.  The data are converted to a
    linear array (could be a list/tuple/array/sequence, I'm not saying),
    *series*, where series[m] gives the temperature (a floating
    point value in degrees C) for month *m*, counting from 0 for the
    January of the first year with data.

    (*series*,*begin*) is returned, where *begin* is
    the first year for which there is data for the station.

    Invalid data are marked in the input file with -9999 but are
    translated in the data arrays to BAD.
    """

    begin = None
    # Year from previous line.
    prev = None
    # The previous line itself.
    prevline = None
    series = []
    for line in lines:
        year = int(line[12:16])
        if begin is None:
            begin = year
        if prev == year:
            # There is one case where there are multiple lines for the
            # same year for a particular station.  The v2.mean input
            # file has 3 identical lines for "8009991400101971"
            if line == prevline:
                print "NOTE: repeated record found: Station %s year %s; data are identical" % (line[:12],line[12:16])
                continue
            # This is unexpected.
            assert 0, "Two lines specify different data for %s" % line[:16]
        # Check that the sequence of years increases.
        assert not prev or prev < year

        # In v2.mean the sequence of years for a station record is not
        # necessarily contiguous.  For example "1486284000001988" is
        # immediately followed by "1486284000001990", missing out 1989.
        # Extend with blank data.
        while prev and prev < year-1:
            series.extend([BAD]*12)
            prev += 1
        prev = year
        prevline = line
        for m in range(12):
            datum = int(line[16+5*m:21+5*m])
            if datum == -9999:
                datum = BAD
            else:
                # Convert to floating point and degrees C.
                datum *= 0.1
            series.append(datum)
    return (series, begin)

def monthly_annual(data):
    """Computing a mean and set of annual anomalies.  This has a
    particular algorithm (via monthly and then seasonal means and
    anomalies), which must be maintained for bit-for-bit compatibility
    with GISTEMP; maybe we can drop it later.
    """

    years = len(data) // 12
    # For each of the 12 months, its mean.
    monthly_mean = array.array('d',[BAD]*12)
    # In month major order, every datum converted to anomaly by
    # subtracting the mean for its month.
    monthly_anom = []
    for m in range(12):
        # *row* contains the data for one month of the year.
        row = data[m::12]
        sum,count = sum_valid(row)
        if count > 0:
            monthly_mean[m] = sum / float(count)
        mean = monthly_mean[m]
        monthly_anom_row = array.array('d',[BAD]*years)
        for n in range(years):
            datum = row[n]
            if valid(mean) and valid(datum):
                monthly_anom_row[n] = datum - mean
        monthly_anom.append(monthly_anom_row)

    # :todo:
    # The seasonal calculation shoud be abstracted into a function.
    # inputs: years, monthly_anom, monthly_mean
    # outputs: seasonal_anom, seasonal_mean
    # Average monthly means to make seasonal means,
    # and monthly anomalies to make seasonal anomalies.
    seasonal_mean = array.array('d',[BAD, BAD, BAD, BAD])
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
        seasonal_anom_row = array.array('d', [BAD] * years)
        for n in range(years):
            sum = 0.0
            count = 0
            for m in months:
                if m == 11:
                    if n == 0:
                        continue
                    else:
                        m_anom = monthly_anom[m][n-1]
                else:
                    m_anom = monthly_anom[m][n]
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
        annual_mean = BAD
    annual_anom = array.array('d',[BAD]*years)
    for n in range(years):
        sum,count = sum_valid(seasonal_anom[s][n] for s in range(4))
        # need 3 valid seasons for a valid year
        if count > 2:
            annual_anom[n] = sum / float(count)
    return (annual_mean, annual_anom)


def read_v2():
    """A generator to iterate through the temperature data in
    'work/v2.mean_comb'.  Each item in the iterator is a (dict,
    series) pair for a single 12-character station ID.
    """

    info = read_config.v2_get_info()
    sources = read_config.v2_get_sources()
    f = open('work/v2.mean_comb', 'r')
    print "reading work/v2.mean_comb"
    def id12(l):
        """Extract the 12-digit station record identifier."""
        return l[:12]

    for (id, lines) in itertools.groupby(f, id12):
        # lines is a set of lines which all begin with the same 12 character id
        series, begin = from_lines(list(lines))
        dict = info[id[:11]].copy()
        dict['begin'] = begin
        dict['source'] = sources.get(id, 'UNKNOWN')
        dict['id' ] = id
        yield (dict, round_series(series))

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

    data = [BAD] * (years*12)
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

    rec_begin = record['dict']['begin']
    rec_years = record['years']
    rec_data = record['data']
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

    new_data = fresh_average(sums, wgts, years)
    ann_mean, ann_anoms = monthly_annual(new_data)
    overlap = 0
    # :todo: the records are consulted in an essentially arbitrary
    # order (chosen by the implementation of items()), but the order
    # may affect the result.
    # Tie breaks go to the last record consulted.
    for rec_id, record in records.items():
        rec_ann_anoms = record['ann_anoms']
        rec_ann_mean = record['ann_mean']
        rec_years = record['years']
        rec_begin = record['dict']['begin']
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
        return 0, 0, BAD
    return best_record, best_id, diff

def combine(sums, wgts, begin, years, records):
    while records:
        record, rec_id, diff = get_longest_overlap(sums, wgts, begin, years, records)
        if invalid(diff):
            comb_log.write("\tno other records okay\n")
            return
        del records[rec_id]
        add(sums, wgts, diff, begin, record)
        rec_begin = record['dict']['begin']
        comb_log.write("\t %s %d %d %f\n" %
          (rec_id, rec_begin, record['years'] + rec_begin - 1, diff))

def get_best(records):
    """Given a set of records (a dict really), return the best one, and
    its key in the *records* dict.
    """

    ranks = {'MCDW': 4, 'USHCN2': 3, 'SUMOFDAY': 2, 'UNKNOWN': 1}
    best = 1
    longest = 0
    for rec_id in sorted(records.keys()):
        record = records[rec_id]
        source = record['dict']['source']
        length = record['length']
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
        (dict, series) = records[rec_id]
        begin = dict['begin']
        years = len(series) // 12
        end = begin + years - 1
        y_min = min(y_min, begin)
        y_max = max(y_max, end)
        ann_mean, ann_anoms = monthly_annual(series)
        # Let *length* be the number of valid data in ann_anoms.
        length = 0
        for anom in ann_anoms:
            length += valid(anom)
        record_dict[rec_id] = {'dict': dict, 'data': series,
                               'years': years, 'length': length,
                               'ann_anoms': ann_anoms, 'ann_mean': ann_mean}
    return record_dict, y_min, y_max

def fresh_arrays(record, begin, years):
    """Make and return a fresh set of sums, and wgts arrays.  Each
    array is list (of length 12 * years).

    *begin* should be the starting year for the arrays, which must
    be no later than the starting year for the record.
    """

    nmonths = years * 12

    rec_data = record['data']
    rec_begin, rec_years = record['dict']['begin'], record['years']
    # Number of months in record.
    rec_months = rec_years * 12
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

def get_id11(item):
    """Given an item from the stream passed to comb_records,
    comb_pieces or similar, return the 11-digit identifier."""

    # Each item is a (meta,series) pair.
    meta,series_ = item
    return meta['id'][:11]

def comb_records(stream):
    """Combine records for the same station (the same id11) where
    possible.  For each station, the number of combined records will be
    at most the number of original records (in the case where no
    combining is possible), and each combined record is yielded."""

    global comb_log
    comb_log = open('log/comb.log','w')

    for id11, record_set in itertools.groupby(stream, get_id11):
        comb_log.write('%s\n' % id11)
        records = {}
        for (dict, series) in record_set:
            records[dict['id']] = (dict, series)
        ids = records.keys()
        while 1:
            if len(ids) == 1:
                yield records[ids[0]]
                break
            record_dict, begin, end = make_record_dict(records, ids)
            years = end - begin + 1
            record, rec_id = get_best(record_dict)
            rec_dict = record['dict']
            del record_dict[rec_id]
            sums, wgts = fresh_arrays(record, begin, years)
            new_dict = rec_dict.copy()
            source = rec_dict['source']
            comb_log.write ("\t%s %s %s -- %s\n" % (rec_id, begin, end,source))
            combine(sums, wgts, begin, years, record_dict)
            begin,final_data = final_average(sums, wgts, years, begin)
            new_dict['begin'] = begin
            yield (new_dict, round_series(final_data))
            ids = record_dict.keys()
            if not ids:
                break

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
        return BAD
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
        rec_ann_anoms = record['ann_anoms']
        rec_years = record['years']
        rec_begin = record['dict']['begin']
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
        # diff = sum / wgt
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
    rec_begin = record['dict']['begin']
    rec_end = rec_begin + record['years'] - 1

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

    rec_ann_anoms = record['ann_anoms']
    rec_ann_mean = record['ann_mean']
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
                    anom1 = BAD
                else:
                    anom1 = new_ann_anoms[index1]
                if index2 < 0 or index2 >= rec_len:
                    anom2 = BAD
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
        rec_begin = record['dict']['begin']
        rec_end = rec_begin + record['years'] - 1

        pieces_log.write("\t %s %d %d\n" % (rec_id, rec_begin, record['years'] + rec_begin - 1))

        is_okay = find_quintuples(sums, wgts,
                                  begin, years, record, rec_begin,
                                  new_id, rec_id)

        if is_okay:
            del records[rec_id]
            add(sums, wgts, 0.0, begin, record)
            pieces_log.write("\t %s %d %d\n" % (rec_id, rec_begin, record['years'] + rec_begin - 1))
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
        length = record['length']
        if length > longest:
            longest = length
            longest_rec = record
            longest_id = rec_id
    return longest_rec, longest_id

def comb_pieces(stream):
    global pieces_log
    pieces_log = open('log/pieces.log','w')
    helena_ds = read_config.get_helena_dict()

    for id11, record_set in itertools.groupby(stream, get_id11):
        pieces_log.write('%s\n' % id11)
        records = {}
        for (dict, series) in record_set:
            records[dict['id']] = (dict, series)
        ids = records.keys()
        while 1:
            if len(ids) == 1:
                yield records[ids[0]]
                break
            record_dict, begin, end = make_record_dict(records, ids)

            if helena_ds.has_key(id11):
                id1, this_year, month, summand = helena_ds[id11]
                dict = record_dict[id1]['dict']
                data = record_dict[id1]['data']
                begin = dict['begin']
                years = record_dict[id1]['years']
                # Index of month specified by helena_ds
                M = (this_year - begin)*12 + month
                # All valid data up to and including M get adjusted
                for i in range(M+1):
                    datum = data[i]
                    if invalid(datum):
                        continue
                    data[i] += summand
                del helena_ds[id11]
                ann_mean, ann_anoms = monthly_annual(data)
                record_dict[id1]['ann_anoms'] = ann_anoms
                record_dict[id1]['ann_mean'] = ann_mean

            record, rec_id = get_longest(record_dict)
            years = end - begin + 1
            rec_dict = record['dict']
            source = rec_dict['source']
            del record_dict[rec_id]
            sums, wgts = fresh_arrays(record, begin, years)
            new_dict = rec_dict.copy()
            pieces_log.write("\t%s %s %s -- %s\n" % (rec_id, begin, begin + years -1,source))
            pieces_combine(sums, wgts, begin, years, record_dict, rec_id)
            begin, final_data = final_average(sums, wgts, years, begin)
            new_dict['begin'] = begin
            yield (new_dict, final_data)
            ids = record_dict.keys()
            if not ids:
                break

def drop_strange(data):
    """Drops data from station records, under control of the file
    'config/Ts.strange.RSU.list.IN' file.  Returns an iterator.
    """

    changes_dict = read_config.get_changes_dict()
    for (dict, series) in data:
        changes = changes_dict.get(dict['id'], [])
        begin = dict['begin']
        end = begin + (len(series)//12) - 1
        for (kind, year, x) in changes:
            if kind == 'years':
                # omit all the data from year1 to year2, inclusive
                year1 = year
                year2 = x
                if year1 <= begin and year2 >= end:
                    # drop this whole record
                    series = None
                    break
                if year1 <= begin:
                    if year2 >= begin:
                        # trim at the start
                        series = series[(year2 - begin + 1)*12:]
                        begin = year2 + 1
                        dict['begin'] = begin
                    continue
                if year2 >= end:
                    if year1 <= end:
                        # trim at the end
                        series = series[:(year1 - begin)*12]
                        end = year1 - 1
                    continue
                # remove some years from mid-series
                nmonths = (year2 + 1 - year1) * 12
                series[(year1-begin)*12:(year2+1-begin)*12] = [BAD]*nmonths
            else: # remove a single month
                series[(year-begin)*12 + x-1] = BAD
        if series is not None:
            yield (dict, series)


def alter_discont(data):
    """Modifies records as specified in config/Ts.discont.RS.alter.IN,
    by adding the delta to every datum for that station prior to the
    specified month.  Returns an iterator.
    """

    alter_dict = read_config.get_alter_dict()
    for (dict, series) in data:
        if alter_dict.has_key(dict['id']):
            (a_month, a_year, a_num) = alter_dict[dict['id']]
            begin = dict['begin']
            # Month index of the month in the config file.
            M = (a_year - begin)*12 + a_month - 1
            # Every (valid) month up to and not including the month in
            # question is adjusted.
            for i in range(M):
                if valid(series[i]):
                    series[i] += a_num
        yield(dict, series)

def to_text(dict, series):
    begin = dict['begin']
    t = ' %4d%5d%s%4s%-36s\n' % (math.floor(dict['lat']*10+0.5),
                                 math.floor(dict['lon']*10+0.5),
                                 dict['id'], dict['elevs'], dict['name'])
    for y in range(len(series)//12):
        l = map(lambda m: "%5d" % math.floor(m * 10.0 + 0.5),
          series[y*12:(y+1)*12])
        t += '%4d%s\n' % (y + begin, ''.join(l))
    return t

def write_to_file(data, file):
    """Writes *data* (an iterable of (dict, series) pairs) to the file
    named *file*.
    """

    f = open(file, 'w')
    for (dict, series) in data:
        f.write(to_text(dict, series))
    f.close()

def results():
    data = read_v2()
    records = comb_records(data)
    combined_pieces = comb_pieces(records)
    without_strange = drop_strange(combined_pieces)
    return alter_discont(without_strange)

def main():
    write_to_file(results(), 'work/Ts.txt')

if __name__ == '__main__':
    main()

# Notes:
#
# 1. read_v2() reads work/v2.mean_comb, produced by step 0.
#
# 2. comb_records() is a filter, combining several separate "scribal"
# station records into a single record, based on a minimum overlap
# (number of years with both stations valid) between pairs of records,
# attempting to make single coherent records.
#
# 3. comb_pieces() is another filter, which further attempts to
# combine the records produced by comb_records() - which therefore
# have shorter overlaps - by comparing the annual anomalies of the
# years in which they do overlap, and finding ones for which the
# temperatures (in years which they do have in common) are on average
# closer together than the standard deviation of the combined record.
# I think.
# 
# Before performing this combination, comb_pieces also makes
# discontinuity adjustments to particular records, under control of a
# configuration file.  In fact, the configuration file specifies a
# single adjustment to a single record (St Helena).
# 
# 4. drop_strange() is another filter, which discards some station
# records, or parts of records, under control of a configuration file.
# 
# 5. alter_discont() makes discontinuity adjustments to records under
# control of a configuration file.  Yes, this is very similar to the
# St Helena behaviour of comb_pieces().
#
# 6. write_to_file() writes the results to work/Ts.txt

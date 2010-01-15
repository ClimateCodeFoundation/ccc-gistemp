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

import math, struct, array

def FEQUALS(x, y):
  return abs(x - y) < 0.1

# For each station record we carry around a dict and a series of
# temperature records.  The series is a 2D array of floating-point
# temperatures: data[m][y], where m ranges from 0 to 11 and y from 0
# to N.  Missing values are recorded as 999.9.
#
# This "2D array" is actually a list of array.array values. For
# instance, data[2] is an array of all the values for March.
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
#
# The rotation between year-major data (in the inputs and outputs
# of this step) and month-major data (everywhere inside this step) is
# done by from_lines() and write_to_file().

def from_lines(lines):
    """*lines* is a list of lines (strings) that comprise a station's
    entire record.  The lines are an extract from a file in the same
    format as the GHCN file v2.mean.Z.  The data are converted to a 2D
    array, *series*, where series[m][yi] gives the temperature (a floating
    point value in degrees C) for month *m* of year begin+yi (*begin* is
    the first year for which there is data for the station).
    (*series*,*begin*) is returned.

    Invalid data are marked in the input file with -9999 but are
    translated in the data arrays to 9999*0.1 (which is not quite
    999.9).
    """

    begin = 9999
    end = 0
    for line in lines:
      year = int(line[12:16])
      if year < begin:
        begin = year
      if year > end:
        end = year
    years = end - begin + 1
    l = [999.9]*years
    series = map(lambda x: array.array('d',l), range(12))
    for line in lines:
      year = int(line[12:16])
      index = year-begin
      for m in range(12):
        datum = int(line[16+5*m:21+5*m])
        if datum == -9999:
          datum = 9999
        series[m][index] = datum * 0.1
    return (series, begin)

# Computing a mean and set of annual anomalies.  This has a particular
# algorithm (via monthly and then seasonal means and anomalies, with
# particular thresholds for missing values at each step), which must
# be maintained for bit-for-bit compatibility with GISTEMP; maybe we
# can drop it later.

def monthly_annual(data):
        years = len(data[0])
        # compute monthly means and anomalies   
        monthly_means = array.array('d',[BAD]*12)
        monthly_anoms = []
        for m in range(12):
            row = data[m]
            sum = 0.0
            count = 0
            for n in range(years):
                datum = row[n]
                if not FEQUALS(datum, BAD):
                    sum += datum
                    count += 1
            if count > 0:
                monthly_means[m] = sum/count
            mean = monthly_means[m]
            monthly_anoms_row = array.array('d',[BAD]*years)
            for n in range(years):
                datum = row[n]
                if (not FEQUALS(mean, BAD) and not FEQUALS(datum, BAD)):
                    monthly_anoms_row[n] = datum - mean
            monthly_anoms.append(monthly_anoms_row)

        seasonal_means = array.array('d',[BAD, BAD, BAD, BAD])
        seasonal_anoms = []
        # compute seasonal anomalies
        for s in range(4):
            if s == 0:
                months = [11, 0, 1]
            else:
                months = range(s*3-1, s*3+2)
            sum = 0.0
            count = 0
            for i in range(3):
                m = months[i]
                monthly_mean = monthly_means[m]
                if not FEQUALS(monthly_mean, BAD):
                    sum += monthly_mean
                    count += 1
            if count > 1:
                seasonal_means[s] = sum / count
            seasonal_anoms_row = array.array('d', [BAD] * years)
            for n in range(years):
                sum = 0.0
                count = 0
                for i in range(3):
                    m = months[i]
                    if m == 11:
                      if n == 0:
                        continue
                      else:
                        monthly_anom = monthly_anoms[m][n-1]
                    else:
                        monthly_anom = monthly_anoms[m][n]
                    if not FEQUALS(monthly_anom, BAD):
                        sum += monthly_anom
                        count += 1
                if count > 1:
                    seasonal_anoms_row[n] = sum / count
            seasonal_anoms.append(seasonal_anoms_row)
        
        # annual mean and anomalies
        sum = 0.0
        count = 0
        for s in range(4):
            seasonal_mean = seasonal_means[s]
            if not FEQUALS(seasonal_mean, BAD):
                sum += seasonal_mean
                count += 1
        if count > 2:
                annual_mean = sum / count
        else:
                annual_mean = BAD
        annual_anoms = array.array('d',[BAD]*years)
        for n in range(years):
            sum = 0.0
            count = 0
            for s in range(4):
                anom = seasonal_anoms[s][n]
                if not FEQUALS(anom, BAD):
                    sum += anom
                    count += 1
            if count > 2:
                annual_anoms[n] = sum / count
        return (annual_mean, annual_anoms)

def v2_get_sources():
    """Reads the three tables mcdw.tbl, ushcn2.tbl, sumofday.tbl and
    return a dictionary that maps from 12-digit (string) station ID to
    the source (which is one of the strings 'MCDW', 'USHCN2',
    'SUMOFDAY').
    """

    sources = {}
    for source in ['MCDW', 'USHCN2', 'SUMOFDAY']:
        for line in open('input/%s.tbl' % source.lower()):
            _, id, rec_no = line.split()
            sources[id + rec_no] = source
    return sources

def v2_get_info():
    """Open the input/v2.inv file and convert it into a dictionary, *d*,
    that is returned.

    d[id] contains the metadata of the station with identification *id*
      where *id* is an 11-digit string.

    The input file is in the same format as the GHCN V2 file
    v2.temperature.inv (in fact, it's the same file, but with records
    added for the Antarctic stations that GHCN doesn't have).  The best
    description of that file's format is the Fortran program:
    ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/v2/v2.read.inv.f
    """

    keys = 'id name lat lon elevs elevg pop ipop topo stveg stloc iloc airstn itowndis grveg'.split()
    info = {}
    for row in open('input/v2.inv', 'r'):
        values = struct.unpack('11sx30sx6sx7sx4s5sc5s2s2s2s2sc2s16sx', row[:101])
        # Would use struct.unpack_from if we knew we had Python 2.5
        hash = dict(zip(keys, values))
        hash['ipop'] = int(hash['ipop'])
        hash['lat'] = float(hash['lat'])
        hash['lon'] = float(hash['lon'])
        info[hash['id']] = hash
    return info

def read_v2():
    """A generator to iterate through the temperature data in
    'work/v2.mean_comb'.  Each item in the iterator is a pair (id11,
    dict), where id11 is the 11-digit station ID and dict is a map
    from 12-digit station ID to a (dict, series) pair.
    """

    f = open('work/v2.mean_comb', 'r')
    print "reading work/v2.mean_comb"
    info = v2_get_info()
    sources = v2_get_sources()
    ids = []
    line = f.readline()
    iterdict = {}
    id11 = line[:11]
    while line:
        lines = [line]
        id = line[:12]
        ids.append(id)
        last_id = id
        line = f.readline()
        while line:
            id = line[:12]
            if id != last_id:
                break
            lines.append(line)
            line = f.readline()
        series, begin = from_lines(lines)
        st_id = last_id[:11]
        dict = info[st_id].copy()
        dict['begin'] = begin
        dict['source'] = sources.get(last_id, 'UNKNOWN')
        dict['id' ] = last_id
        if last_id[:11] != id11:
          yield (id11, iterdict)
          iterdict = {}
          id11 = last_id[:11]
        iterdict[last_id] = (dict, series)
    yield (id11, iterdict)

BAD = 9999 / 10.0
BAD = 0.1 * int(BAD * 10.0)
MIN_OVERLAP = 4
comb_log = None

def average(new_sums, new_counts, new_data, years):
    for m in range(12):
        for n in range(years):
            count = new_counts[m][n]
            if count != 0:
                new_data[m][n] = new_sums[m][n] / count

def final_average(new_sums, new_wgts, new_data, years, begin):
    min, max = 9999, -9999
    for m in range(12):
        mon_min, mon_max = 9999, -9999
        wgts_row = new_wgts[m]
        for n in range(years):
            wgt = wgts_row[n]
            if wgt == 0:
                continue
            if mon_min > n:
                mon_min = n
            if mon_max < n:
                mon_max = n
        if min > mon_min:
            min = mon_min
        if max < mon_max:
            max = mon_max
    if min == 0 and max == years - 1:
        average(new_sums, new_wgts, new_data, years)
        return begin
    years = max - min + 1
    begin = begin + min
    end = begin + years - 1
    for m in range(12):
        new_sums[m] = new_sums[m][min: max + 1]
        new_wgts[m] = new_wgts[m][min: max + 1]
        new_data[m] = [BAD] * years
    average(new_sums, new_wgts, new_data, years)
    return begin

def add(new_sums, new_wgts, diff, begin, record):
    rec_begin = record['dict']['begin']
    rec_years = record['years']
    rec_data = record['data']
    for m in range(12):
        sums_row = new_sums[m]
        wgts_row = new_wgts[m]
        data_row = rec_data[m]
        for n in range(rec_years):
            datum = data_row[n]
            if abs(datum - BAD) < 0.1:
                continue
            year = n + rec_begin
            index = year - begin
            sums_row[index] = sums_row[index] + datum - diff
            wgts_row[index] = wgts_row[index] + 1

def get_longest_overlap(new_sums, new_wgts, new_data, begin, years, records):
    average(new_sums, new_wgts, new_data, years)
    ann_mean, ann_anoms = monthly_annual(new_data)
    overlap = 0
    length = 0
    for rec_id, record in records.items():
        rec_ann_anoms = record['ann_anoms']
        rec_ann_mean = record['ann_mean']
        rec_years = record['years']
        rec_begin = record['dict']['begin']
        sum = wgt = 0
        for n in range(rec_years):
            rec_anom = rec_ann_anoms[n]
            if abs(rec_anom - BAD) < 0.1:
                continue
            year = n + rec_begin
            anom = ann_anoms[year - begin]
            if abs(anom - BAD) < 0.1:
                continue
            wgt = wgt + 1
            sum = sum + (rec_ann_mean + rec_anom) - (ann_mean + anom)
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

def combine(new_sums, new_wgts, new_data, begin, years, records):
    while records:
        record, rec_id, diff = get_longest_overlap(new_sums, new_wgts, new_data, begin, years, records)
        if abs(diff - BAD) < 0.1:
            comb_log.write("\tno other records okay\n")
            return
        del records[rec_id]
        add(new_sums, new_wgts, diff, begin, record)
        rec_begin = record['dict']['begin']
        comb_log.write("\t %s %d %d %f\n" % (rec_id, rec_begin, record['years'] + rec_begin - 1, diff))

def get_best(records):
    ranks = {'MCDW': 4, 'USHCN2': 3, 'SUMOFDAY': 2, 'UNKNOWN': 1}
    best = 1
    longest = 0
    rids = records.keys()
    rids.sort()
    for rec_id in rids:
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
    record_dict = {}
    min, max = 9999, -9999
    for rec_id in ids:
        s = records[rec_id]
        (dict, data) = s
        begin = dict['begin']
        years = len(data[0])
        end = begin + years - 1
        if min > begin:
            min = begin
        if max < end:
            max = end
        ann_mean, ann_anoms = monthly_annual(data)
        length = 0
        for anom in ann_anoms:
            if anom != BAD:
                length = length + 1
        record_dict[rec_id] = {'dict': dict, 'data': data,
                               'string': s, 'years': years,
                               'length': length, 'ann_anoms': ann_anoms,
                               'ann_mean': ann_mean}
    return record_dict, min, max

def records_get_new_data(record, begin, years):
    sums, wgts, data = [None] * 12, [None] * 12, [None] *12
    rec_data = record['data']
    rec_begin, rec_years = record['dict']['begin'], record['years']
    rec_end = rec_begin + rec_years - 1
    for m in range(12):
        sums_row, wgts_row, data[m] = [0.0] * years, [0] * years, [BAD] * years
        rec_row = rec_data[m]
        for n in range(rec_years):
            datum = rec_row[n]
            if abs(datum - BAD) < 0.1:
                continue
            index = n + rec_begin - begin
            sums_row[index] = datum
            wgts_row[index] = 1
        sums[m] = sums_row
        wgts[m] = wgts_row
    return sums, wgts, data

def comb_records(stream):
    global comb_log
    comb_log = open('log/comb.log','w')

    for (id, records) in stream:
        comb_log.write('%s\n' % id)
        iterdict = {}
        ids = records.keys()
        while 1:
            record_dict, begin, end = make_record_dict(records, ids)
            if len(record_dict) == 1:
                rec_id, record = record_dict.items()[0]
                iterdict[rec_id] = (record['dict'],record['data'])
                yield (id, iterdict)
                break
            years = end - begin + 1
            record, rec_id = get_best(record_dict)
            rec_dict = record['dict']
            source = rec_dict['source']
            del record_dict[rec_id]
            new_sums, new_wgts, new_data = records_get_new_data(record, begin, years)
            new_dict = rec_dict.copy()
            comb_log.write ("\t%s %s %s -- %s\n" % (rec_id, begin, end,source))
            combine(new_sums, new_wgts, new_data, begin, years, record_dict)
            begin = final_average(new_sums, new_wgts, new_data, years, begin)
            new_dict['begin'] = begin
            iterdict[rec_id] = (new_dict, new_data)
            ids = record_dict.keys()
            if not ids:
              yield (id, iterdict)
              break
    comb_log.close()

MIN_MID_YEARS = 5            # MIN_MID_YEARS years closest to ymid
BUCKET_RADIUS = 10
pieces_log = None

def sigma(list, bad):
    sos = sum = count = 0
    for x in list:
        if x == bad: continue
        sos = sos + x * x
        sum = sum + x
        count = count + 1
    if count == 0:
        return bad
    else:
        mean = sum/count
        return math.sqrt(sos / count - mean * mean)

def pieces_get_longest_overlap(new_sums, new_wgts, new_data, begin, years, records):
    average(new_sums, new_wgts, new_data, years)
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
            if abs(rec_anom - BAD) < 0.1:
                continue
            year = n + rec_begin
            anom = ann_anoms[year - begin]
            if abs(anom - BAD) < 0.1:
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

def get_actual_endpoints(new_wgts, begin, years):
    tmp_begin = 9999
    tmp_end = 0
    for m in range(12):
        wgts_row = new_wgts[m]
        for n in range(years):
            if wgts_row[n] == 0:
                continue
            year = n + begin
            if year < tmp_begin:
                tmp_begin = year
            if tmp_end < year:
                tmp_end = year
    return tmp_begin, tmp_end

def find_quintuples(new_sums, new_wgts, new_data,
                    begin, years, record, rec_begin,
                    new_id, rec_id):
    rec_begin = record['dict']['begin']
    rec_end = rec_begin + record['years'] - 1

    actual_begin, actual_end = get_actual_endpoints(new_wgts, begin, years)

    max_begin = max(actual_begin, rec_begin)
    min_end = min(actual_end, rec_end)
    middle_year = int(.5 * (max_begin + min_end) + 0.5)
    pieces_log.write("max begin: %s\tmin end: %s\n" % (max_begin, min_end))

    average(new_sums, new_wgts, new_data, years)
    ann_mean, sublist1 = monthly_annual(new_data)
    ann_std_dev = sigma(sublist1, BAD)
    pieces_log.write("ann_std_dev = %s\n" % ann_std_dev)
    offset1 = (middle_year - begin)
    len1 = len(sublist1)

    sublist2 = record['ann_anoms']
    rec_ann_mean = record['ann_mean']
    offset2 = (middle_year - rec_begin)
    len2 = len(sublist2)

    ov_success = 0
    okay_flag = 0
    for rad in range(1, BUCKET_RADIUS + 1):
        count1 = sum1 = 0
        count2 = sum2 = 0
        for i in range(0, rad + 1):
            for sign in [-1, 1]:
                if sign == 1 and i == 0:
                    continue
                index1 = i * sign + offset1
                index2 = i * sign + offset2
                if index1 < 0 or index1 >= len1:
                    anom1 = BAD
                else:
                    anom1 = sublist1[index1]
                if index2 < 0 or index2 >= len2:
                    anom2 = BAD
                else:
                    anom2 = sublist2[index2]
                if anom1 != BAD:
                    sum1 = sum1 + anom1 + ann_mean
                    count1 = count1 + 1
                if anom2 != BAD:
                    sum2 = sum2 + anom2 + rec_ann_mean
                    count2 = count2 + 1
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

def pieces_combine(new_sums, new_wgts, new_data, begin, years, records, new_id):
    while records:
        record, rec_id = pieces_get_longest_overlap(new_sums, new_wgts, new_data, begin, years, records)
        rec_begin = record['dict']['begin']
        rec_end = rec_begin + record['years'] - 1

        pieces_log.write("\t %s %d %d\n" % (rec_id, rec_begin, record['years'] + rec_begin - 1))

        is_okay = find_quintuples(new_sums, new_wgts, new_data,
                                  begin, years, record, rec_begin,
                                  new_id, rec_id)

        if is_okay:
            del records[rec_id]
            add(new_sums, new_wgts, 0.0, begin, record)
            pieces_log.write("\t %s %d %d\n" % (rec_id, rec_begin, record['years'] + rec_begin - 1))
        else:
            pieces_log.write("\t***no other pieces okay***\n")
            return

def get_longest(records):
    longest = 0
    for rec_id, record in records.items():
        length = record['length']
        if length > longest:
            longest = length
            longest_rec = record
            longest_id = rec_id
    return longest_rec, longest_id

def get_helena_dict():
    """Reads the file config/combine_pieces_helena.in into a dict,
    mapping a station id to a tuple (ID with duplicate marker, year,
    month, summand)."""
    
    helena_ds = {}
    for line in open('config/combine_pieces_helena.in', 'r'):
        id, _, year, month, summand = line.split()
        helena_ds[id[:-1]] = (id, int(year), int(month), float(summand))
    return helena_ds

def comb_pieces(stream):
    global pieces_log
    pieces_log = open('log/pieces.log','w')
    helena_ds = get_helena_dict()

    for (id, records) in stream:
        pieces_log.write('%s\n' % id)
        ids = records.keys()
        while 1:
            record_dict, begin, end = make_record_dict(records, ids)

            if len(record_dict) == 1:
                rec_id, record = record_dict.items()[0]
                yield (record['dict'], record['data'])
                break

            if helena_ds.has_key(id):
                id1, this_year, month, summand = helena_ds[id]
                dict = record_dict[id1]['dict']
                data = record_dict[id1]['data']
                begin = dict['begin']
                years = record_dict[id1]['years']
                for n in range(years):
                    year = n + begin
                    if year > this_year:
                        break
                    for m in range(12):
                        if year == this_year and m > month:
                            break
                        datum = data[m][n]
                        if datum == BAD:
                            continue
                        data[m][n] = datum + summand
                del helena_ds[id]
                ann_mean, ann_anoms = monthly_annual(data)
                record_dict[id1]['ann_anoms'] = ann_anoms
                record_dict[id1]['ann_mean'] = ann_mean

            record, rec_id = get_longest(record_dict)
            years = end - begin + 1
            rec_dict = record['dict']
            source = rec_dict['source']
            del record_dict[rec_id]
            new_sums, new_wgts, new_data = records_get_new_data(record, begin, years)
            new_dict = rec_dict.copy()
            pieces_log.write("\t%s %s %s -- %s\n" % (rec_id, begin, begin + years -1,source))
            pieces_combine(new_sums, new_wgts, new_data, begin, years, record_dict, rec_id)
            begin = final_average(new_sums, new_wgts, new_data, years, begin)
            new_dict['begin'] = begin
            yield (new_dict, new_data)
            ids = record_dict.keys()
            if not ids:
                break

def get_changes_dict():
    """Reads the file config/Ts.strange.RSU.list.IN and returns a dict
    result.  Each line in that file begins with a 12-digit station ID
    - actually the tuple (country-code, WMO station, modifier,
    duplicate) - and ends with either yyyy/mm, specifying a month
    datum to omit or with xxxx-yyyy, specifying years to omit.  xxxx
    can be 0, meaning from the beginning. yyyy can be 9999, meaning to
    the end.  The dict is a map from ID to ('month',yyyy,mm) or
    ('years',xxxx,yyyy).
    """

    dict = {}
    for line in open('config/Ts.strange.RSU.list.IN', 'r'):
        split_line = line.split()
        id = split_line[0]
        try:
            year1, year2 = map(int, split_line[-1].split("-"))
            val = ("years", year1, year2)
        except ValueError:
            year, month = map(int, split_line[-1].split("/"))
            val = ("month", year, month)
        dict[id] = dict.get(id,[])
        dict[id].append(val)
    return dict

def drop_strange(data):
    """Drops data from station records, under control of the file
    'config/Ts.strange.RSU.list.IN' file.  A filter generator.
    """

    changes_dict = get_changes_dict()
    for (dict, series) in data:
        changes = changes_dict.get(dict['id'], [])
        begin = dict['begin']
        end = begin + len(series[0]) - 1
        for (kind, year, x) in changes:
            if kind == 'years': # omit all the data from year1 to year2, inclusive
                year1 = year
                year2 = x
                if (year1 <= begin and year2 >= end): # drop this whole record
                    series = None
                    break
                if (year1 <= begin):
                    if (year2 >= begin): # trim at the start
                        for m in range(12):
                            series[m] = series[m][year2 - begin + 1:]
                        begin = year2 + 1
                        dict['begin'] = begin
                    continue
                if (year2 >= end):
                    if (year1 <= end): # trim at the end
                        for m in range(12):
                            series[m] = series[m][:year1 - begin]
                        end = year1 - 1
                    continue
                # remove some years from mid-series
                for year in range(year1, year2 + 1):
                    for m in range(12):
                        series[m][year - begin] = BAD
            else: # remove a single month
                series[x-1][year-begin] = BAD
        if series is not None:
            yield (dict, series)

def get_alter_dict():
    """Reads the file config/Ts.discont.RS.alter.IN into a dict.  Each
    line has a 12 digit station ID, a month, a year, and a
    floating-point temperature delta.  The dict maps the ID to (month,
    year, delta).
    """

    dict = {}
    for line in open('config/Ts.discont.RS.alter.IN'):
        id, month, year, num = line.split()
        dict[id] = [int(month), int(year), float(num)]
    return dict

def alter_discont(data):
    """Modifies records as specified in config/Ts.discont.RS.alter.IN,
    by adding the delta to every datum for that station prior to the
    specified month.  A filter generator.
    """

    alter_dict = get_alter_dict()
    for (dict, series) in data:
        if (alter_dict.has_key(dict['id'])):
            (a_month, a_year, a_num) = alter_dict[dict['id']]
            begin = dict['begin']
            for year in range(begin, a_year + 1):
                index = year - begin
                for m in range(12):
                    if year == a_year and m > a_month - 2:
                        continue
                    if series[m][index] != BAD:
                        series[m][index] += a_num
        yield(dict, series)

def to_text(dict, series):
    begin = dict['begin']
    t = ' %4d%5d%s%4s%-36s\n' % (math.floor(dict['lat']*10+0.5),
                                 math.floor(dict['lon']*10+0.5),
                                 dict['id'], dict['elevs'], dict['name'])
    for y in range(len(series[0])):
        l = map(lambda m: "%5d" % math.floor(series[m][y] * 10.0 + 0.5), range(12))
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
# 1. read_v2() reads work/v2.mean_comb
#
# 2. comb_records() is a filter
#
# 3. comb_pieces() is a filter
#
# 4. drop_strange() is a filter
#
# 5. alter_discont() is a filter
#
# 6. write_to_file() writes the results to work/Ts.txt
